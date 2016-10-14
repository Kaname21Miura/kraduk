# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 16:53:36 2016

@author: Kaname
"""

#装置のパーツを格納したモジュール
import math
import random

# スリットの親クラス ####
class Slit(object):
    def __init__(self,outerD,slitD,width,thickness,position):
        self.outerD = outerD
        self.slitD = slitD
        self.width = width
        self.thickness = thickness
        self.position = position
        self.d_out = self.slitD/2 + self.width/2
        self.d_in = self.slitD/2 - self.width/2
        self.position = position
        self.front_z = position + thickness
        self.back_z = position
    
    def hittingPotision(self,x,y,z,ux,uy,uz,hit_z):
        db = (hit_z - z)/uz
        x_slit = x + db*ux
        y_slit = y + db*uy
        z_slit = hit_z
        return x_slit,y_slit,z_slit


#### レンズの親クラス ####
class Lens (object):
    def __init__(self,outerD,efl,bfl,ct,et,r,n,position,typ):
        self.outerD = outerD
        self.efl = efl
        self.bfl = bfl
        self.ct = ct
        self.et = et
        self.r = r
        self.n = n
        self.position = position
        if typ == 1:
            self.center = self.position - (self.ct - self.r)
        else:
            if typ == 2:
                self.center = self.position + (self.ct - self.r)
            else:
                print("レンズの向きが入力されていないか間違っています")
        self.x_l = 0
        self.y_l = 0
        self.z_l = self.center

    def hittingPointPlano(self,x,y,z,ux,uy,uz,hit_z):
        db = (hit_z - z)/uz
        x = x + db*ux
        y = y + db*uy
        z = hit_z
        return x,y,z

    def hittingPointConvex(self,x,y,z,ux,uy,uz,typ): #円弧の内側からの衝突の場合 typ = 1,外側からの場合 typ = 2
        delta = ((x-self.x_l)*ux+(y-self.y_l)*uy+(z-self.z_l)*uz)**2-(x-self.x_l)**2-(y-self.y_l)**2-(z-self.z_l)**2+self.r**2
        if typ == 1:
            t = -((x-self.x_l)*ux + (y-self.y_l)*uy + (z-self.z_l)*uz) + math.sqrt(delta)
        if typ == 2:
            t = -((x-self.x_l)*ux + (y-self.y_l)*uy + (z-self.z_l)*uz) - math.sqrt(delta)
        x = x + t*ux
        y = y + t*uy
        z = z + t*uz
        return x,y,z

    def VectorConvPlano(self,x,y,z,ux,uy,uz,e,direction):
        nx,ny,nz = self.orthVectorPlano()
        snel_co = self.snelLow(direction)
        ux = snel_co*(ux - (ux*nx + uy*ny + uz*nz)*nx) + nx*math.sqrt(e)
        uy = snel_co*(uy - (ux*nx + uy*ny + uz*nz)*ny) + ny*math.sqrt(e)
        uz = snel_co*(uz - (ux*nx + uy*ny + uz*nz)*nz) + nz*math.sqrt(e)
        return ux,uy,uz

    def VectorConvConvex(self,x,y,z,ux,uy,uz,e,typ,direction):
        nx,ny,nz = self.orthVectorConvex(x,y,z,typ)
        snel_co = self.snelLow(direction)
        ux = snel_co*(ux - (ux*nx + uy*ny + uz*nz)*nx) + nx*math.sqrt(e)
        uy = snel_co*(uy - (ux*nx + uy*ny + uz*nz)*ny) + ny*math.sqrt(e)
        uz = snel_co*(uz - (ux*nx + uy*ny + uz*nz)*nz) + nz*math.sqrt(e)
        return ux,uy,uz

    def reflect(self,ux,uy,uz,nx,ny,nz):
        e = 1 - ((1/self.n)**2)*(1 - (ux*nx + uy*ny + uz*nz))
        return e

    def orthVectorPlano(self):
        nx = 0
        ny = 0
        nz = -1
        return nx,ny,nz

    def orthVectorConvex(self,x,y,z,typ):
        if typ == 1: #球の内側から
            nx = (x-self.x_l)/self.r
            ny = (y-self.y_l)/self.r
            nz = (z-self.z_l)/self.r
        if typ == 2: #球の外側から
            nx = -(x-self.x_l)/self.r
            ny = -(y-self.y_l)/self.r
            nz = -(z-self.z_l)/self.r
        return nx,ny,nz

    def snelLow(self,direction):
        if direction == 1:
            snel_co = 1/self.n
        if direction == 0:
            snel_co = self.n/1
        return snel_co

    def reflectance(self,ux,uy,uz,old_ux,old_uy,old_uz,nx,ny,nz,nt): #内部反射率
        ai = self.angle(old_ux,old_uy,old_uz,nx,ny,nz)
        if (self.n > nt) and (ai > math.asin(nt/self.n)) :
            Ra = 1
        else:
            at = self.angle(ux,uy,uz,nx,ny,nz)
            rs = (math.sin(ai-at)/math.sin(ai+at))**2
            rp = (math.tan(ai-at)/math.tan(ai+at))**2
            Ra = (rs+rp)/2
        return Ra

    def angle(self,ux,uy,uz,nx,ny,nz): #入射角,屈折角
        ai = math.acos((ux*nx + uy*ny + uz*nz)/(math.sqrt(ux**2 + uy**2 + uz**2)*math.sqrt(nx**2 + ny**2 + nz**2)))
        return ai



##### 曲面が負の方向を向いているレンズ #####
##内部反射率導入モデル##
class Lens1(Lens):
    def __init__(self,outerD,efl,bfl,ct,et,r,n,position):
        super().__init__(outerD,efl,bfl,ct,et,r,n,position,1)
    
    def opticalAnalysis(self,x,y,z,ux,uy,uz):
        x,y,z = self.hittingPointPlano(x,y,z,ux,uy,uz,self.position)
        nx,ny,nz = self.orthVectorPlano()
        e = self.reflect(ux,uy,uz,nx,ny,nz)
        if e <= 0:
            leReflect = True
        else:
            old_ux = ux; old_uy = uy; old_uz = uz
            ux,uy,uz = self.VectorConvPlano(x,y,z,ux,uy,uz,e,1)
            Ra = self.reflectance(ux,uy,uz,old_ux,old_uy,old_uz,nx,ny,nz,self.n)
            if random.random() <= Ra:
                leReflect = True
            else:
                
                x,y,z = self.hittingPointConvex(x,y,z,ux,uy,uz,1)
                nx,ny,nz = self.orthVectorConvex(x,y,z,1)
                e = self.reflect(ux,uy,uz,nx,ny,nz)
                if e <= 0:
                    leReflect = True
                else:
                    old_ux = ux; old_uy = uy; old_uz = uz
                    ux,uy,uz = self.VectorConvConvex(x,y,z,ux,uy,uz,e,1,0)
                    Ra = self.reflectance(ux,uy,uz,old_ux,old_uy,old_uz,nx,ny,nz,1)
                    if random.random() <= Ra:
                        leReflect = True
                    else:
                        leReflect = False
                    
        return x,y,z,ux,uy,uz,leReflect

##### 曲面が正の方向を向いたレンズ #####
##内部反射率導入モデル##
class Lens2(Lens):
    def __init__(self,outerD,efl,bfl,ct,et,r,n,position):
        super().__init__(outerD,efl,bfl,ct,et,r,n,position,2)
    
    def opticalAnalysis(self,x,y,z,ux,uy,uz):
        x,y,z = self.hittingPointConvex(x,y,z,ux,uy,uz,2)
        nx,ny,nz = self.orthVectorConvex(x,y,z,2)
        e = self.reflect(ux,uy,uz,nx,ny,nz)
        if e <= 0:
            leReflect = True
        else:
            old_ux = ux; old_uy = uy; old_uz = uz
            ux,uy,uz = self.VectorConvConvex(x,y,z,ux,uy,uz,e,2,1)
            Ra = self.reflectance(ux,uy,uz,old_ux,old_uy,old_uz,nx,ny,nz,self.n)
            if random.random() <= Ra:
                leReflect = True
            else:
                
                x,y,z = self.hittingPointPlano(x,y,z,ux,uy,uz,self.position)
                nx,ny,nz = self.orthVectorPlano()
                e = self.reflect(ux,uy,uz,nx,ny,nz)
                if e <= 0:
                    leReflect = True
                else:
                    old_ux = ux; old_uy = uy; old_uz = uz
                    ux,uy,uz = self.VectorConvPlano(x,y,z,ux,uy,uz,e,0)
                    Ra = self.reflectance(ux,uy,uz,old_ux,old_uy,old_uz,nx,ny,nz,1)
                    if random.random() <= Ra:
                        leReflect = True
                    else:
                        leReflect = False
                    
        return x,y,z,ux,uy,uz,leReflect


#### フォトダイオードのクラス ####
class Photodiode(object):
    def __init__(self,d,position):
        self.d = d
        self.r = d/2
        self.position = position
        self.count = 0
        self.record_w = 0
    
    def hittingPotision(self,x,y,z,ux,uy,uz,hit_z):
        db = (hit_z - z)/uz
        x = x + db*ux
        y = y + db*uy
        z = hit_z
        return x,y,z
    
    def absorbance(self,nPh):
        if self.record_w == 0:
            ab = 0
        else:
            ab = math.log(self.record_w/nPh)
        return ab



"""
#### レンズの親クラス ####
class Lens (object):
    def __init__(self,outerD,efl,bfl,ct,et,r,n,position,typ):
        self.outerD = outerD
        self.efl = efl
        self.bfl = bfl
        self.ct = ct
        self.et = et
        self.r = r
        self.n = n
        self.position = position
        if typ == 1:
            self.center = self.position - (self.ct - self.r)
        else:
            if typ == 2:
                self.center = self.position + (self.ct - self.r)
            else:
                print("レンズの向きが入力されていないか間違っています")
        self.x_l = 0
        self.y_l = 0
        self.z_l = self.center

    def hittingPointPlano(self,x,y,z,ux,uy,uz,hit_z):
        db = (hit_z - z)/uz
        x = x + db*ux
        y = y + db*uy
        z = hit_z
        return x,y,z

    def hittingPointConvex(self,x,y,z,ux,uy,uz,typ): #円弧の内側からの衝突の場合 typ = 1,外側からの場合 typ = 2
        delta = ((x-self.x_l)*ux+(y-self.y_l)*uy+(z-self.z_l)*uz)**2-(x-self.x_l)**2-(y-self.y_l)**2-(z-self.z_l)**2+self.r**2
        if typ == 1:
            t = -((x-self.x_l)*ux + (y-self.y_l)*uy + (z-self.z_l)*uz) + math.sqrt(delta)
        if typ == 2:
            t = -((x-self.x_l)*ux + (y-self.y_l)*uy + (z-self.z_l)*uz) - math.sqrt(delta)
        x = x + t*ux
        y = y + t*uy
        z = z + t*uz
        return x,y,z

    def VectorConvPlano(self,x,y,z,ux,uy,uz,e,direction):
        nx,ny,nz = self.orthVectorPlano()
        snel_co = self.snelLow(direction)
        ux = snel_co*(ux - (ux*nx + uy*ny + uz*nz)*nx) + nx*math.sqrt(e)
        uy = snel_co*(uy - (ux*nx + uy*ny + uz*nz)*ny) + ny*math.sqrt(e)
        uz = snel_co*(uz - (ux*nx + uy*ny + uz*nz)*nz) + nz*math.sqrt(e)
        return ux,uy,uz

    def VectorConvConvex(self,x,y,z,ux,uy,uz,e,typ,direction):
        nx,ny,nz = self.orthVectorConvex(x,y,z,typ)
        snel_co = self.snelLow(direction)
        ux = snel_co*(ux - (ux*nx + uy*ny + uz*nz)*nx) + nx*math.sqrt(e)
        uy = snel_co*(uy - (ux*nx + uy*ny + uz*nz)*ny) + ny*math.sqrt(e)
        uz = snel_co*(uz - (ux*nx + uy*ny + uz*nz)*nz) + nz*math.sqrt(e)
        return ux,uy,uz

    def reflect(self,ux,uy,uz,nx,ny,nz):
        e = 1 - ((1/self.n)**2)*(1 - (ux*nx + uy*ny + uz*nz))
        return e

    def orthVectorPlano(self):
        nx = 0
        ny = 0
        nz = -1
        return nx,ny,nz

    def orthVectorConvex(self,x,y,z,typ):
        if typ == 1:
            nx = (x-self.x_l)/self.r
            ny = (y-self.y_l)/self.r
            nz = (z-self.z_l)/self.r
        if typ == 2:
            nx = -(x-self.x_l)/self.r
            ny = -(y-self.y_l)/self.r
            nz = -(z-self.z_l)/self.r
        return nx,ny,nz

    def snelLow(self,direction):
        if direction == 1:
            snel_co = 1/self.n
        if direction == 0:
            snel_co = self.n/1
        return snel_co


##### 曲面が負の方向を向いているレンズ #####
class Lens1(Lens):
    def __init__(self,outerD,efl,bfl,ct,et,r,n,position):
        super().__init__(outerD,efl,bfl,ct,et,r,n,position,1)
    
    def opticalAnalysis(self,x,y,z,ux,uy,uz):
        x,y,z = self.hittingPointPlano(x,y,z,ux,uy,uz,self.position)
        nx,ny,nz = self.orthVectorPlano()
        e = self.reflect(ux,uy,uz,nx,ny,nz)
        if e <= 0:
            leReflect = True
        else:
            leReflect = False
            ux,uy,uz = self.VectorConvPlano(x,y,z,ux,uy,uz,e,1)
            
            x,y,z = self.hittingPointConvex(x,y,z,ux,uy,uz,1)
            nx,ny,nz = self.orthVectorConvex(x,y,z,1)
            e = self.reflect(ux,uy,uz,nx,ny,nz)
            if e <= 0:
                leReflect = True
            else:
                leReflect = False
                ux,uy,uz = self.VectorConvConvex(x,y,z,ux,uy,uz,e,1,0)
        return x,y,z,ux,uy,uz,leReflect

##### 曲面が正の方向を向いたレンズ #####
class Lens2(Lens):
    def __init__(self,outerD,efl,bfl,ct,et,r,n,position):
        super().__init__(outerD,efl,bfl,ct,et,r,n,position,2)
    
    def opticalAnalysis(self,x,y,z,ux,uy,uz):
        x,y,z = self.hittingPointConvex(x,y,z,ux,uy,uz,2)
        nx,ny,nz = self.orthVectorConvex(x,y,z,2)
        e = self.reflect(ux,uy,uz,nx,ny,nz)
        if e <= 0:
            leReflect = True
        else:
            leReflect = False
            ux,uy,uz = self.VectorConvConvex(x,y,z,ux,uy,uz,e,2,1)
            
            x,y,z = self.hittingPointPlano(x,y,z,ux,uy,uz,self.position)
            nx,ny,nz = self.orthVectorPlano()
            e = self.reflect(ux,uy,uz,nx,ny,nz)
            if e <= 0:
                leReflect = True
            else:
                leReflect = False
                ux,uy,uz = self.VectorConvPlano(x,y,z,ux,uy,uz,e,0)
        return x,y,z,ux,uy,uz,leReflect
"""

