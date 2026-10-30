import cadquery as cq, math

# ---------------- parameters (mm) ----------------
R_DISC=125.0; DISC_T=10.0
R_HAT_O=67.0; R_HAT_I=59.5
PLATE_FRONT=-14.0; PLATE_T=3.0          # key plate (flange of MX switch rests on its front)
FLANGE_PCBZ=6.5                          # switch flange height above PCB bottom (from STEP)
dz=PLATE_FRONT-FLANGE_PCBZ               # world z = pcb z + dz  (board bottom)
BOARD_BOT=dz; BOARD_TOP=dz+1.51
HAT_BACK=-36.0; COVER_Z0,COVER_Z1=-36.0,-33.0
KEY_HOLE=14.2
SCR_W,SCR_H,SCR_T=62.5,17.9,2.3; ACT_W,ACT_H=55.29,14.8
SCR_Y=15.0
BOARD_CX,BOARD_CY=155.5,-101.5
SW_X=[-29.04,-10.04,9.0,28.08]; SW_Y=-8.0
MOUNT=[(44,21),(-44,21),(44,-21),(-44,-21)]
SCREW_ANG=[30,150,270]

pcb=cq.importers.importStep('/mnt/user-data/uploads/PCB.step')
pcb_w=cq.Workplane("XY").add(pcb.val().translate((-BOARD_CX,-BOARD_CY,dz)))

def cyl(r,h,x,y,z): return cq.Workplane("XY").add(cq.Solid.makeCylinder(r,h,cq.Vector(x,y,z),cq.Vector(0,0,1)))
def box(w,d,h,x,y,z): return cq.Workplane("XY").box(w,d,h).translate((x,y,z))
def rad(r,h,ang,z,rin):                  # radial cylinder starting at radius rin
    d=cq.Vector(math.cos(math.radians(ang)),math.sin(math.radians(ang)),0)
    return cq.Workplane("XY").add(cq.Solid.makeCylinder(r,h,cq.Vector(rin*d.x,rin*d.y,z),d))

# ---------------- rotor front: disc + hat + key plate ----------------
ring=cyl(R_DISC,DISC_T,0,0,-DISC_T).cut(cyl(R_HAT_O-1,DISC_T,0,0,-DISC_T))
wall=cyl(R_HAT_O,-HAT_BACK,0,0,HAT_BACK).cut(cyl(R_HAT_I,-HAT_BACK,0,0,HAT_BACK))
plate=cyl(R_HAT_I+0.1,PLATE_T,0,0,PLATE_FRONT-PLATE_T)
rotor=ring.union(wall).union(plate)
rotor=rotor.cut(cyl(74,1.0,0,0,-1.0).cut(cyl(72,1.0,0,0,-1.0)))      # decorative step groove

# screen pad (front of plate) with pocket, window, cable slot
PAD_T=3.0
pad=box(SCR_W+3.6,SCR_H+3.6,PAD_T,0,SCR_Y,PLATE_FRONT+PAD_T/2)
rotor=rotor.union(pad)
pocket_z_top=PLATE_FRONT+PAD_T
rotor=rotor.cut(box(SCR_W+0.4,SCR_H+0.4,SCR_T+0.2,0,SCR_Y,pocket_z_top-(SCR_T+0.2)/2+0.1))
rotor=rotor.cut(box(ACT_W+1.0,ACT_H+1.0,PAD_T+PLATE_T+2,0,SCR_Y,PLATE_FRONT-PLATE_T/2+PAD_T/2))
rotor=rotor.cut(box(3.0,15.4,PAD_T+PLATE_T+2,-(SCR_W+0.4)/2-1.0,SCR_Y,PLATE_FRONT-PLATE_T/2+PAD_T/2))
for x in SW_X: rotor=rotor.cut(box(KEY_HOLE,KEY_HOLE,PLATE_T+2,x,SW_Y,PLATE_FRONT-PLATE_T/2))

# curved vent holes
pts=[]
for g in range(5):
    for rowoff in (-3.2,3.2):
        for i in range(11):
            r=79+i*3.8
            th=math.radians(g*72+15+rowoff*100/r+0.9*(r-79))
            pts.append((r*math.cos(th),r*math.sin(th)))
for (x,y) in pts: rotor=rotor.cut(cyl(2.1,DISC_T+1,x,y,-DISC_T-0.5))

# USB-C panel slot (+Y) at XIAO USB height, and radial cover screws
USB_Z=BOARD_BOT-2.15
rotor=rotor.cut(box(12,12,6,2.45,63,USB_Z))
for a in SCREW_ANG: rotor=rotor.cut(rad(1.7,14,a,-34.5,55))

# ---------------- back cover with standoffs ----------------
cover=cyl(R_HAT_I-0.25,COVER_Z1-COVER_Z0,0,0,COVER_Z0)
for (x,y) in MOUNT:
    cover=cover.union(cyl(3.5,BOARD_BOT-COVER_Z1,x,y,COVER_Z1)).cut(cyl(1.3,10,x,y,BOARD_BOT-10))
for a in SCREW_ANG: cover=cover.cut(rad(1.3,6,a,-34.5,R_HAT_I-0.25).translate((0,0,0)).mirror("XY") if False else
    cq.Workplane("XY").add(cq.Solid.makeCylinder(1.3,6,cq.Vector((R_HAT_I-0.25)*math.cos(math.radians(a)),(R_HAT_I-0.25)*math.sin(math.radians(a)),-34.5),cq.Vector(-math.cos(math.radians(a)),-math.sin(math.radians(a)),0))))

# ---------------- stand ----------------
stand=box(150,35,52,0,-122.5,-20).cut(cyl(R_DISC+0.2,12,0,0,-11))

def vol(a,b):
    try: return round(a.intersect(b).val().Volume(),2)
    except Exception: return 0.0
print("PCB/rotor overlap:",vol(pcb_w,rotor)," PCB/cover:",vol(pcb_w,cover)," rotor/stand:",vol(rotor,stand)," cover/rotor:",vol(cover,rotor))
bb=pcb_w.val().BoundingBox(); print("PCB z range",round(bb.zmin,1),round(bb.zmax,1))

out="/mnt/user-data/outputs/"
for n,o in (("rotor_front",rotor),("back_cover",cover),("stand",stand)):
    cq.exporters.export(o,out+n+".step"); cq.exporters.export(o,out+n+".stl",tolerance=0.1)
asm=cq.Assembly()
asm.add(rotor,name="rotor_front",color=cq.Color(0.6,0.6,0.62)); asm.add(cover,name="back_cover",color=cq.Color(0.3,0.3,0.32))
asm.add(stand,name="stand",color=cq.Color(0.25,0.15,0.1)); asm.add(pcb_w,name="pcb",color=cq.Color(0.1,0.5,0.3))
asm.export(out+"rotor_clock_assembly.step")
print("done")
