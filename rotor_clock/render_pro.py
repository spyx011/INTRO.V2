import numpy as np, numba as nb, pickle, trimesh, math, sys
from scipy.ndimage import uniform_filter, gaussian_filter
from PIL import Image

P=pickle.load(open('/home/claude/parts.pkl','rb'))
# material: (r,g,b, metallic, clearcoat/spec strength, F0)
MAT={"rotor":(0.74,0.75,0.78,0.92,1.0,0.60),"cal":(0.96,0.20,0.015,0.0,1.0,0.06),"bolt":(0.90,0.90,0.93,1.0,1.0,0.85),
     "cover":(0.13,0.13,0.15,0.35,0.9,0.12),"usb":(0.80,0.81,0.84,1.0,1.0,0.7),"tft":(0.02,0.02,0.025,0.0,0.8,0.05),
     "glass":(0.003,0.004,0.007,0.0,1.0,0.14),"pcb":(0.05,0.33,0.15,0.0,0.6,0.04),"sw":(0.05,0.05,0.06,0.0,0.9,0.06),"cable":(0.03,0.03,0.035,0.0,0.7,0.05),"carrier":(0.20,0.20,0.23,0.2,0.8,0.1)}
NAMES=list(MAT.keys()); MATARR=np.array([MAT[k] for k in NAMES],dtype=np.float32)

def prep_mesh(V,F):
    m=trimesh.Trimesh(V,F,process=False); m.merge_vertices()
    try:
        sm=m.smooth_shaded
        return np.asarray(sm.vertices),np.asarray(sm.faces),np.asarray(sm.vertex_normals)
    except Exception:
        return np.asarray(m.vertices),np.asarray(m.faces),np.asarray(m.vertex_normals)
CACHE={}
def get(name):
    if name not in CACHE:
        V,F=P[name]; CACHE[name]=prep_mesh(V,F)
    return CACHE[name]

@nb.njit(cache=True)
def raster(X,Y,Z,N,F,mid,zbuf,nbuf,mbuf):
    H,W=zbuf.shape
    for t in range(F.shape[0]):
        a=F[t,0];b=F[t,1];c=F[t,2]
        x0=X[a];y0=Y[a];x1=X[b];y1=Y[b];x2=X[c];y2=Y[c]
        area=(x1-x0)*(y2-y0)-(x2-x0)*(y1-y0)
        if abs(area)<1e-9: continue
        minx=max(0,int(math.floor(min(x0,min(x1,x2)))));maxx=min(W-1,int(math.ceil(max(x0,max(x1,x2)))))
        miny=max(0,int(math.floor(min(y0,min(y1,y2)))));maxy=min(H-1,int(math.ceil(max(y0,max(y1,y2)))))
        for py in range(miny,maxy+1):
            cy=py+0.5
            for px in range(minx,maxx+1):
                cx=px+0.5
                w0=((x1-cx)*(y2-cy)-(x2-cx)*(y1-cy))/area
                w1=((x2-cx)*(y0-cy)-(x0-cx)*(y2-cy))/area
                w2=1.0-w0-w1
                if w0<-1e-5 or w1<-1e-5 or w2<-1e-5: continue
                z=w0*Z[a]+w1*Z[b]+w2*Z[c]
                if z>zbuf[py,px]:
                    zbuf[py,px]=z
                    nbuf[py,px,0]=w0*N[a,0]+w1*N[b,0]+w2*N[c,0]
                    nbuf[py,px,1]=w0*N[a,1]+w1*N[b,1]+w2*N[c,1]
                    nbuf[py,px,2]=w0*N[a,2]+w1*N[b,2]+w2*N[c,2]
                    mbuf[py,px]=mid

def sstep(a,b,x): t=np.clip((x-a)/(b-a),0,1); return t*t*(3-2*t)
LIGHTS=[(np.array([-0.55,0.65,0.55]),5.0,0.20,np.array([1.0,0.98,0.95])),
        (np.array([0.80,0.20,0.45]),2.6,0.26,np.array([0.85,0.92,1.0])),
        (np.array([0.0,-0.5,0.6]),0.8,0.35,np.array([1.0,1.0,1.0])),
        (np.array([0.0,1.0,0.25]),2.0,0.30,np.array([1.0,1.0,1.0])),
        (np.array([0.3,0.45,-0.85]),2.2,0.25,np.array([1.0,1.0,1.0]))]
for i,(d,_,_,_) in enumerate(LIGHTS): LIGHTS[i]=(d/np.linalg.norm(d),)+LIGHTS[i][1:]
def env(d):
    up=d[...,1:2]
    sky=(0.50+0.0*up)*np.array([1.0,1.03,1.08])
    zen=np.array([0.95,1.0,1.12]); gnd=np.array([0.04,0.04,0.05]); hor=np.array([0.30,0.32,0.37])
    top=hor+(zen-hor)*sstep(0,0.9,up); bot=hor+(gnd-hor)*sstep(0,0.7,-up)
    c=np.where(up>=0,top,bot)
    for (L,I,w,col) in LIGHTS:
        dp=np.sum(d*L,axis=-1,keepdims=True); c=c+I*col*np.exp(-(1-dp)/(w*w*0.5))
    return c

def render(fn,az,el,hide=(),lim=150,cx=0,cy=6,W=1100,ss=2,seed=0):
    a=math.radians(az); e=math.radians(el)
    Ry=np.array([[math.cos(a),0,math.sin(a)],[0,1,0],[-math.sin(a),0,math.cos(a)]])
    Rx=np.array([[1,0,0],[0,math.cos(e),-math.sin(e)],[0,math.cos(e)*0+math.sin(e)*0+1,0]])  # placeholder (overwritten below)
    Rx=np.array([[1,0,0],[0,math.cos(e),-math.sin(e)],[0,math.sin(e),math.cos(e)]])
    M=Rx@Ry
    Wp=W*ss; Hp=W*ss; sc=Wp/(2*lim)
    zbuf=np.full((Hp,Wp),-1e9,np.float32); nbuf=np.zeros((Hp,Wp,3),np.float32); mbuf=np.full((Hp,Wp),-1,np.int8)
    for name in P:
        if name in hide: continue
        V,F,Nn=get(name)
        v=V@M.T; n=Nn@M.T
        X=((v[:,0]-cx)*sc+Wp/2).astype(np.float64); Y=(Hp/2-(v[:,1]-cy)*sc).astype(np.float64); Z=v[:,2].astype(np.float64)
        raster(X,Y,Z,n.astype(np.float64),F.astype(np.int64),NAMES.index(name),zbuf,nbuf,mbuf)
    hit=mbuf>=0
    # screen-space ambient occlusion from the depth buffer
    zb=np.where(hit,zbuf,np.nan); zmin=np.nanmin(zb); zf=np.where(hit,zbuf,zmin)
    k=int(0.045*Wp)|1
    blur=uniform_filter(zf,size=k)
    wgt=uniform_filter(hit.astype(np.float32),size=k)+1e-6
    blur=uniform_filter(np.where(hit,zbuf,0),size=k)/wgt
    ao=np.clip((blur-zf)/(10*sc/ (Wp/(2*lim))*0+10.0),0,1)   # depth difference in mm
    ao=1-0.60*np.clip(ao*1.1,0,1)**0.8
    out=np.zeros((Hp,Wp,3),np.float32)
    rows=range(0,Hp,300)
    for r0 in rows:
        r1=min(Hp,r0+300); m=mbuf[r0:r1]; N=nbuf[r0:r1].copy()
        nl=np.linalg.norm(N,axis=-1,keepdims=True)+1e-9; N=N/nl
        flip=N[...,2:3]<0; N=np.where(flip,-N,N)
        V_=np.array([0,0,1],np.float32); NV=np.clip(N[...,2:3],0,1)
        R=2*NV*N-V_
        mi=np.clip(m,0,None); mp=MATARR[mi]
        alb=mp[...,0:3]; metal=mp[...,3:4]; spec=mp[...,4:5]; F0=mp[...,5:6]
        # diffuse: hemispherical ambient + lights
        amb=env(N*0.0+np.array([0,1,0]))*0.0
        hemi=0.20+0.16*N[...,1:2]
        dif=hemi*np.array([1.0,1.0,1.04])
        for (L,I,w,col) in LIGHTS[:3]:
            dif=dif+0.16*I*col*np.clip(np.sum(N*L,axis=-1,keepdims=True),0,1)
        fres=F0+(1-F0)*(1-NV)**5
        E=env(R)
        diffcol=alb*dif*(1-metal)
        speccol=E*(alb*metal+(1-metal)*1.0)*fres*spec*np.where(metal>0.5,1.0,1.0)
        metalbase=alb*metal*E*0.88
        col=diffcol+speccol*(1-metal)+metalbase+metal*speccol*0.25
        col=col*ao[r0:r1,:,None]
        out[r0:r1]=col
    # filmic tone map + gamma
    out=out/(1+0.40*out); out=np.clip(out,0,1)**(1/2.2)
    cov=hit.astype(np.float32)
    out=out*cov[...,None]
    # downsample
    o=out.reshape(W,ss,W,ss,3).mean(axis=(1,3)); c=cov.reshape(W,ss,W,ss).mean(axis=(1,3))
    yy,xx=np.mgrid[0:W,0:W]; rr=np.hypot(xx-W/2,yy-W*0.46)/(W*0.75)
    bg=(0.965-0.20*np.clip(rr,0,1)**1.6)[...,None]*np.array([0.98,0.99,1.0])
    img=o+(1-c[...,None])*bg
    Image.fromarray((np.clip(img,0,1)*255).astype(np.uint8)).save(fn)

if __name__=="__main__":
    view=sys.argv[1] if len(sys.argv)>1 else "all"
    if view in("all","iso"):   render('/home/claude/iso.png',-33,10,hide=())
    if view in("all","front"): render('/home/claude/front.png',0,0,hide=())
    if view in("all","back"):  render('/home/claude/back.png',152,12,hide=("pcb",))
    if view in("all","hub"):   render('/home/claude/hub.png',-18,6,hide=(),lim=82,cy=0)
    if view in("all","inside"): render('/home/claude/inside.png',-22,28,hide=("cover","tft","glass","rotor","cal"),lim=75,cy=8)
    if view in("all","cal"):   render('/home/claude/caliper.png',-8,-6,hide=(),lim=78,cx=0,cy=98)
