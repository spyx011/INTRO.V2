import cadquery as cq, math, numpy as np, trimesh, tempfile, pickle, os

STEP='/mnt/user-data/uploads/BREAK_DISC_Assembly.STEP'
PCBSTEP='/mnt/user-data/uploads/PCB.step'
OUT='/mnt/user-data/outputs/rotor_clock/'
K=1.25                                   # uniform scale so the PCB fits the hub (rotor geometry otherwise untouched)
CY,CZ,X0=871.5,1364.84,724.8             # disc axis (y,z) and hat inner-floor plane (x) in the source file

# ---------------- place the ORIGINAL rotor into the clock frame ----------------
# clock frame: Z = toward the viewer (looking into the hub), Y = up (caliper at the bottom), origin on the hub's inner floor
asm=cq.importers.importStep(STEP); sols=asm.solids().vals()
def place(sol):
    s=sol.translate(cq.Vector(-X0,-CY,-CZ)).scale(K)
    s=s.rotate(cq.Vector(0,0,0),cq.Vector(0,1,0),-90)
    return s.rotate(cq.Vector(0,0,0),cq.Vector(0,0,1),180)     # caliper at the BOTTOM
disc,hat,cal=[place(s) for s in (sols[0],sols[1],sols[2])]
for n,s in (("disc",disc),("hat",hat),("caliper",cal)):
    b=s.BoundingBox(); print("%-8s x %.1f..%.1f  y %.1f..%.1f  z %.2f..%.2f"%(n,b.xmin,b.xmax,b.ymin,b.ymax,b.zmin,b.zmax))
Z_FLOOR_OUT=hat.BoundingBox().zmin          # back of hat (-4.375)
Z_LEDGE=disc.BoundingBox().zmin             # hat rim / disc front plane (14.375)
Z_REAR=disc.BoundingBox().zmax              # disc rear plane = front of clock (30.0)
R_CUP_IN=46.5*K; R_BORE=50.0*K; R_LUG=24.0*K
print("cup inner R %.2f  disc bore R %.2f  ledge z %.2f  rear z %.2f"%(R_CUP_IN,R_BORE,Z_LEDGE,Z_REAR))

# ---------------- PCB / TFT parameters ----------------
CARRIER_T=2.9
Z_BB=CARRIER_T+9.62+1.0                  # board bottom: clears buzzer / header underneath
PLATE_BACK=Z_BB+3.5; PLATE_FRONT=Z_BB+6.5; PLATE_T=3.0     # MX flange rests on the plate front
KEY_HOLE=14.2
SCR_W,SCR_H,SCR_T=62.5,17.9,2.3; ACT_W,ACT_H=55.29,14.8
SCR_Y=15.0
BOARD_CX,BOARD_CY=155.5,-101.5
SW_X=[-29.04,-10.04,9.0,28.08]; SW_Y=-8.0
MOUNT=[(44,21),(-44,21),(44,-21),(-44,-21)]
print("plate front %.2f  key tops %.2f  clock front plane %.2f  recess %.2f"%(PLATE_FRONT,Z_BB+16.8,Z_REAR,Z_REAR-PLATE_FRONT))

pcb=cq.importers.importStep(PCBSTEP)
pcb_w=cq.Workplane("XY").add(pcb.val().translate(cq.Vector(-BOARD_CX,-BOARD_CY,Z_BB)))

def cyl(r,h,x,y,z): return cq.Workplane("XY").add(cq.Solid.makeCylinder(r,h,cq.Vector(x,y,z),cq.Vector(0,0,1)))
def box(w,d,h,x,y,z): return cq.Workplane("XY").box(w,d,h).translate((x,y,z))

# ---------------- carrier insert (sits on the hub floor, holds the PCB, with a plain cable hole) ----------------
carrier=cyl(R_CUP_IN-0.2,CARRIER_T,0,0,0)
lugs=[(R_LUG*math.cos(math.radians(270+72*i)),R_LUG*math.sin(math.radians(270+72*i))) for i in range(5)]
carrier=carrier.cut(cq.Workplane("XY").workplane(offset=-1).pushPoints(lugs).circle(4.5).extrude(CARRIER_T+2))   # keep lug holes open
for (x,y) in MOUNT:
    carrier=carrier.union(cyl(3.5,Z_BB-CARRIER_T,x,y,CARRIER_T)).cut(cyl(1.3,10,x,y,Z_BB-10))
# plain round cable hole at the back: same size as the hat's centre bore (R13 x scale)
HOLE_R=13.0*K
carrier=carrier.cut(cyl(HOLE_R,CARRIER_T+2,0,0,-1))

# ---------------- face plate (keys + screen), rests on the hub ledge ----------------
R_PL=R_BORE-0.2
plate=cyl(R_PL,PLATE_T,0,0,PLATE_BACK)
skirt=cyl(R_PL,PLATE_BACK-Z_LEDGE+0.01,0,0,Z_LEDGE).cut(cyl(R_CUP_IN-1.2,PLATE_BACK-Z_LEDGE+0.02,0,0,Z_LEDGE-0.005))
plate=plate.union(skirt)
try: plate=plate.faces(">Z").edges().chamfer(0.8)
except Exception: pass
PAD_T=3.0
plate=plate.union(box(SCR_W+3.6,SCR_H+3.6,PAD_T,0,SCR_Y,PLATE_FRONT+PAD_T/2).edges("|Z").fillet(3.0).faces(">Z").edges().chamfer(0.6))
top=PLATE_FRONT+PAD_T
plate=plate.cut(box(SCR_W+0.4,SCR_H+0.4,SCR_T+0.2,0,SCR_Y,top-(SCR_T+0.2)/2+0.1))
plate=plate.cut(box(ACT_W+1.0,ACT_H+1.0,PAD_T+PLATE_T+2,0,SCR_Y,PLATE_BACK+(PLATE_T+PAD_T)/2).edges("|Z").fillet(1.5))
plate=plate.cut(box(3.0,15.4,PAD_T+PLATE_T+2,-(SCR_W+0.4)/2-1.0,SCR_Y,PLATE_BACK+(PLATE_T+PAD_T)/2))
for x in SW_X: plate=plate.cut(box(KEY_HOLE,KEY_HOLE,PLATE_T+2,x,SW_Y,PLATE_BACK+PLATE_T/2).edges("|Z").fillet(1.0))
# fine machined grooves for a refined gauge look
for (ra,rb) in ((53.2,54.0),(56.4,57.2)):
    plate=plate.cut(cyl(rb,0.6,0,0,PLATE_FRONT-0.6).cut(cyl(ra,0.6,0,0,PLATE_FRONT-0.6)))
# screw bosses clamp the PCB onto the carrier posts
for (x,y) in MOUNT:
    plate=plate.union(cyl(3.0,PLATE_BACK-(Z_BB+1.51)+0.01,x,y,Z_BB+1.51))
    plate=plate.cut(cyl(1.75,PLATE_BACK-Z_BB+PLATE_T+2,x,y,Z_BB-0.5))
    plate=plate.cut(cq.Workplane("XY").add(cq.Solid.makeCone(1.75,3.2,1.45,cq.Vector(x,y,PLATE_FRONT-1.45),cq.Vector(0,0,1))))

# ---------------- TFT module reference ----------------
tft_body=box(SCR_W,SCR_H,SCR_T,0,SCR_Y,top-SCR_T/2).edges("|Z").fillet(0.6)
tft_glass=box(ACT_W,ACT_H,0.3,0,SCR_Y,top-0.13)

# ---------------- checks ----------------
def cqvol(a,b):
    try: return round(a.intersect(b).val().Volume() if hasattr(a,'val') else a.intersect(b).Volume(),2)
    except Exception as e: return "err"
print("PCB vs carrier",cqvol(pcb_w,carrier)," PCB vs plate",cqvol(pcb_w,plate)," TFT vs plate",cqvol(tft_body,plate))
print("carrier vs hat",cqvol(carrier,cq.Workplane("XY").add(hat))," plate vs hat",cqvol(plate,cq.Workplane("XY").add(hat))," plate vs disc",cqvol(plate,cq.Workplane("XY").add(disc)))
pb=pcb_w.val().BoundingBox(); r=0
for v in pcb_w.val().Vertices(): r=max(r,math.hypot(v.X,v.Y))
print("PCB max radius %.1f vs cup inner %.1f ; PCB z %.1f..%.1f ; key tops %.1f"%(r,R_CUP_IN,pb.zmin,pb.zmax,Z_BB+16.8))

# ---------------- export ----------------
os.makedirs(OUT,exist_ok=True)
def stl(shape,name,tol=0.05):
    w=cq.Workplane("XY").add(shape) if not isinstance(shape,cq.Workplane) else shape
    cq.exporters.export(w,OUT+name,tolerance=tol,angularTolerance=0.1)
for n,s in (("rotor_disc",disc),("rotor_hat",hat),("caliper",cal)):
    cq.exporters.export(cq.Workplane("XY").add(s),OUT+n+".step"); stl(s,n+".stl",0.06)
for n,s in (("carrier",carrier),("face_plate",plate)):
    cq.exporters.export(s,OUT+n+".step"); stl(s,n+".stl",0.02)
cq.exporters.export(tft_body,OUT+"tft_2.25in_reference.step"); pass

def tm(shape,tol=0.06):
    f=tempfile.mktemp(suffix=".stl"); w=cq.Workplane("XY").add(shape) if not isinstance(shape,cq.Workplane) else shape
    cq.exporters.export(w,f,tolerance=tol,angularTolerance=0.1); return trimesh.load(f,process=True)
parts={"rotor":trimesh.util.concatenate([tm(disc),tm(hat,0.03)]),"cal":tm(cal,0.05),"cover":tm(plate,0.02),"carrier":tm(carrier,0.03),
       "tft":tm(tft_body,0.02),"glass":tm(tft_glass,0.02)}
_ps=pcb_w.solids().vals()
parts["sw"]=trimesh.util.concatenate([tm(x,0.05) for x in _ps[:4]])
parts["pcb"]=trimesh.util.concatenate([tm(x,0.05) for x in _ps[4:]])
pickle.dump({k:(np.asarray(v.vertices),np.asarray(v.faces)) for k,v in parts.items()},open('/home/claude/parts.pkl','wb'))
sc=trimesh.Scene(); cols={"rotor":[170,172,180,255],"cal":[235,95,20,255],"cover":[40,40,45,255],"carrier":[45,45,52,255],"tft":[15,15,18,255],"glass":[5,5,8,255],"pcb":[25,120,70,255],"cable":[22,22,26,255],"sw":[20,20,24,255],"usb":[210,210,215,255]}
for k,v in parts.items():
    t=v.copy(); t.visual.face_colors=cols[k]; sc.add_geometry(t,geom_name=k)
sc.export(OUT+"rotor_clock_assembly.glb")
print("exported")

# ---------------- combined STEP assembly (colours) ----------------
A=cq.Assembly(name="rotor_clock")
def add(shape,name,rgb):
    w=shape if isinstance(shape,cq.Workplane) else cq.Workplane("XY").add(shape)
    A.add(w,name=name,color=cq.Color(*rgb))
add(disc,"rotor_disc",(0.66,0.67,0.70)); add(hat,"rotor_hat",(0.66,0.67,0.70)); add(cal,"caliper",(0.92,0.37,0.08))
add(carrier,"carrier",(0.16,0.16,0.18)); add(plate,"face_plate",(0.18,0.18,0.20))
add(pcb_w,"pcb",(0.10,0.45,0.25)); add(tft_body,"tft",(0.05,0.05,0.06)); pass
A.export(OUT+"rotor_clock_assembly.step")
print("assembly step exported")
