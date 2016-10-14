# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 10:57:43 2016

@author: Kaname
"""

#光学的モンテカルロ法のモジュール
#乱数を用いた計算を繰り返すこと近似的に物質中の光子の伝搬挙動を求めていきます

import numpy as np
import math
import random
import time

class Photon(object):
    def __init__ (self):
        self.w = 1                        #光子のエネルギー
        self._wMin = 0.0001               #光子の最小エネルギー
        self.u = np.array([0,0,1])        #光子の方向ベクトル
        self.coor = np.array([0,0,0])  #光子の位置座標
        self.nPh = 0


class Tissue(object): #組織の親クラス
    def __init__(self,g,ms,ma,n,d,z_0):
        self.g = g
        self.ms = ms      #scattering coefficient
        self.ma = ma      #absorption coefficient
        self.mt = ma + ms #interaction coefficient
        self.n = n
        self.d = d
        self.z_0 = z_0
        self.z_1 = z_0 + d
    
    def photonMoving(self,x,y,z,ux,uy,uz,s): #光子の移動
        l = s/self.mt
        nx = x + ux*l
        ny = y + uy*l
        nz = z + uz*l
        return nx,ny,nz
    
    def photonAbsorption(self,w): #光子の吸収
        nw = w * (1 - self.ma/self.mt)
        return nw
    
    def photonScattering(self): #光子の散乱
        #r = random.random()
        shita = math.acos((1/(2*self.g))*(1+self.g**2-((1-(self.g)**2)/(1-self.g+2*self.g*random.random()))**2))
        fai = 2*math.pi*random.random()
        return shita, fai
    
    def vectorConv(self,ux,uy,uz): #ベクトルの変更
        s,f = self.photonScattering()
        nux = math.sin(s)*(ux*uz*math.cos(f) - uy*math.sin(f))/(math.sqrt(1-uz**2)) + ux*math.cos(s)
        nuy = math.sin(s)*(uy*uz*math.cos(f) + ux*math.sin(f))/(math.sqrt(1-uz**2)) + uy*math.cos(s)
        nuz = -math.sin(s)*math.cos(f)*math.sqrt(1-uz**2) + uz*math.cos(s)
        return nux, nuy, nuz
    
    def exVectorConv(self,ux,uy,uz): #|uz| > 0.99999の時に発動
        s,f = self.photonScattering()
        nux = math.sin(s)*math.cos(f)
        nuy = math.sin(s)*math.sin(f)
        nuz = uz*math.cos(s)/abs(uz)
        return nux, nuy, nuz
    
    
    def distanceBoundary(self,uz,z): #境界までの距離
        if uz < 0:
            db = (self.z_0 - z)/uz
        else:
            if uz > 0:
                db = (self.z_1 - z)/uz
            else:
                db = 100
        return db
    
    def hittngPotision(self,x,y,z,ux,uy,uz,db): #境界との衝突地点
        x = x + ux*db
        y = y + uy*db
        z = z + uz*db
        return x, y, z
    
    def angleOfIncidence(self,uz): #入射角
        ai = math.acos(abs(uz))
        return ai
    
    def angleOfTrancemission(self,uz,nt): #屈折角
        ai = self.angleOfIncidence(uz)
        at = math.asin((self.n/nt)*math.sin(ai))
        return at
    
    def reflectance(self,uz,nt): #内部反射率
        ai = self.angleOfIncidence(uz)
        if (self.n > nt) and (ai > math.asin(nt/self.n)) :
            Ra = 1
        else:
            if abs(uz) > 0.9999:
                Ra = 0
            else:
                at = self.angleOfTrancemission(uz,nt)
                rs = (math.sin(ai-at)/math.sin(ai+at))**2
                rp = (math.tan(ai-at)/math.tan(ai+at))**2
                Ra = (rs+rp)/2
        return Ra
    
    def newDirectionByRef(self,ux,uy,uz):
        uz_n = -uz
        return ux,uy,uz_n
    
    def newDirectionByTra(self,ux,uy,uz,nt):
        at = self.angleOfTrancemission(uz,nt)
        ux_n = ux*self.n/nt
        uy_n = uy*self.n/nt
        uz_n = uz*math.cos(at)/abs(uz)
        return ux_n,uy_n,uz_n
        
        
def save(fname,list_sample):
    f = open(fname,'w')
    try:
        for i in range(len(list_sample)):
            for j in range (len(list_sample[i])-1):
                sub_list = round(list_sample[i][j],6)
                f.write(str(sub_list))
                f.write(',')
            sub_list = round(list_sample[i][len(list_sample[i])-1],6)
            f.write(str(sub_list))
            f.write('\n')
            if list_sample[i][4] == 0:
                break
        time.sleep(1.0)
    finally:
        f.close
        

def monte(nPh,fname):
    timer_start = time.time()
    #nPh = 10000000 #光子数
    #fname = '161011Ms600Skin0.0nPh10^7Surf.csv'
    
    
    number = nPh / 10
    
    ###パラメーターの設定###
    #それぞれの大きさはcmで記載
    
    #皮膚
    g_skin = 0.9   #異方性係数
    ms_skin = 150  #散乱係数
    ma_skin = 0.5  #吸収係数
    n_skin = 1.4   #屈折率
    d_skin = 0   #厚さ
    
    #骨
    g_bone = 0.9
    ms_bone = 600
    ma_bone = 0.2
    n_bone = 1.4
    d_bone = 60
    
    #それぞれの位置関係
    z_skin = 0
    z_bone = d_skin
    
    
    skin = Tissue(g_skin,ms_skin,ma_skin,n_skin,d_skin,z_skin)
    bone = Tissue(g_bone,ms_bone,ma_bone,n_bone,d_bone,z_bone)
    
    list_sample = np.zeros((nPh, 7)).astype(np.float32)
    flag_1 = False
    count = 0
    rap_time = 0
    n_count = 0
    rd = 0
    
    while count < 10:
        
        
        sub_number = 0
        
        while number > sub_number:
            sub_number = sub_number + 1
            ph = Photon()
            Rsp = ((1-1.4)**2)/((1+1.4)**2)
            ph.w = ph.w - Rsp
            ux = 0
            uy = 0
            uz = 1
            x = 0
            y = 0
            z = 0
            z_bottom = 0
            
            while 1:
                
                s = -math.log(random.random())
                tissue = bone
                
                if (z >= bone.z_1 ):
                    print("out of tissue")
                    flag_1 = True
                    break
                
                if(d_skin == 0):
                    pass
                else:
                    if (z < skin.d):
                        tissue = skin
                    if (z == skin.d) and (uz < 0):
                        tissue = skin
                
                db = tissue.distanceBoundary(uz, z)
                
                if db * tissue.mt <= s:
                    s = s - db*tissue.mt
                    x,y,z = tissue.hittngPotision(x,y,z,ux,uy,uz,db)
                    
                    if d_skin == 0:
                        nt = 1
                    else:
                        if tissue == bone:
                            nt = skin.n
                        else:
                            if uz > 0:
                                nt = bone.n
                            else:
                                nt = 1
                
                    Ra = tissue.reflectance(uz,nt)
                    
                    if random.random() <= Ra:
                        
                        ux,uy,uz = tissue.newDirectionByRef(ux,uy,uz)
                    else:
                        ux,uy,uz = tissue.newDirectionByTra(ux,uy,uz,nt)
                        if(z <= 0):
                            flag_1 = False
                            rd = rd + ph.w
                            break
    
                else:
        
                    x,y,z = tissue.photonMoving(x,y,z,ux,uy,uz,s)
            
                    if z <= 0:#例外処理
                        db = abs((z-tissue.z_0)/uz)
                        x,y,z = tissue.hittngPotision(x,y,z,ux,uy,uz,db)
                            
                        nt = 1
                        Ra = tissue.reflectance(uz,nt)
                        if random.random() <= Ra:
                                
                            ux,uy,uz = tissue.newDirectionByRef(ux,uy,uz)
                            continue
                        else:
                            ux,uy,uz =tissue.newDirectionByTra(ux,uy,uz,nt)
                                
                            flag_1 = False
                            break
            
                
                    if z_bottom < z:
                        z_bottom = z
                        
                        s = 0
                        
                        ph.w = tissue.photonAbsorption(ph.w)
                        
                        if abs(uz) <= 0.99999:
                            ux,uy,uz = tissue.vectorConv(ux,uy,uz)
                        
                        else:
                            ux,uy,uz = tissue.exVectorConv(ux,uy,uz)
                        
                        ud = math.sqrt(ux**2 + uy**2 + uz**2)
                        ux = ux/ud; uy = uy/ud; uz = uz/ud
    
                    if ph.w <= ph._wMin:
                        flag_1 = True
                        break
                    if flag_1:
                        continue
    
            list_sample[n_count] = [x,y,ux,uy,uz,ph.w,z_bottom]
            n_count = n_count + 1
        
        
        
        rap_time = time.time() - timer_start
        print (rap_time)
        print(count)
        count = count + 1
    
    save(fname,list_sample)
    
    elapsed_time = time.time() - timer_start
    print(elapsed_time)
    print("END")

    


#このモジュールファイルのテストコード
"""if __name__ == '__main__':
    nPn = 10
    fname = 'testMonte.csv'    
    main(nPn,fname)
"""
    
