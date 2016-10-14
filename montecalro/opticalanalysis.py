# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 15:56:43 2016

@author: Kaname
"""

#光学計算モジュール
#光式骨密度計における光路解析を行います　
#その後強度分布を計測します

import numpy as np
import matplotlib.pyplot as plt
import math
import time
import pandas as pa
from os import chdir
from montecalro import opticalparts as parts



def densitometer():
    
    #平凸レンズでの光路計算
    timer_start = time.time()
    workspace = '/Users/Kaname/Documents/ipython'
    #workspaceまでパスを通す
    chdir(workspace)
    
    nPh = 100000 #光子数
    #計算する試料表面のファイル名
    fname =  '160628Ms300Skin0.1nPh10^5surf.csv'
    
    #捜査範囲（ここだけ ｍｍ 記述になります）
    start = -10
    end = 60
    split = 10  #何ｍｍごとに分割するかを入力
    
    step = int((end-start)/split) #総ステップ数
    count = 0
    
    #パラメータの設定(cm)
    #レンズ１
    outerD_1 = 5    #外径
    efl_1 = 10      #焦点距離
    bfl_1 = 9.341   #後焦点距離
    ct_1 = 1        #中心厚
    et_1 = 0.3553   #コバ厚
    r_1 = 5.168     #曲率半径
    n_1 = 1.517     #屈折率
    
    #レンズ２
    outerD_2 = 5
    efl_2 = 5
    bfl_2 = 4.328
    ct_2 = 1.2
    et_2 = 0.301
    r_2 = 3.924
    n_2 = 1.758
    
    #スリット
    slit_outerD = 5  #外径
    slitD = 2        #スリット径
    width = 0.2    #スリット間隙の大きさ
    thickness = 0.3  #厚さ
    
    #フォトダイオード
    d_pd = 0.6       #受光面直径
    
    
    
    
    #その他の変数等
    list_absorbance = np.zeros(step+1)
    list_count = np.zeros(step+1)
    list_w =np.zeros(step+1)
    list_distance = []
    list_zb = []
    list_x = []
    list_y = []
    #Zの位置をlistとして設定
    for i in range (step+1):
        list_distance = list_distance + [start + split*i]
    
    leReflect = False #レンズ境界面での反射(Trueだと反射、Falseだと透過)
    
    #Zのステップごとのループ
    while count <= step:
        out_count = 0
        
        #それぞれの位置関係
        z_record = (start + split*count)*0.1  #z_recordは一回の移動距離を表す
        z_lens1 = -bfl_1 + z_record           #レンズは平面部分を基準に位置を規定
        z_lens2 = z_lens1 - ct_1 -5.5 - ct_2  #レンズ１からコバ厚と６ｃｍ分離した
        z_slit1 = z_lens1                     #スリットの位置はそれぞれレンズの位置と負方向の面を基準として規定
        z_slit2 = z_lens2 + ct_2
        z_pd = z_lens2 - bfl_2
        
        #レンズとスリット、フォトダイオードのオブジェクトをそれぞれ生成
        lens_1 = parts.Lens1(outerD_1,efl_1,bfl_1,ct_1,et_1,r_1,n_1,z_lens1)
        lens_2 = parts.Lens2(outerD_2,efl_2,bfl_2,ct_2,et_2,r_2,n_2,z_lens2)
        slit_1 = parts.Slit(slit_outerD,slitD,width,thickness,z_slit1)
        slit_2 = parts.Slit(slit_outerD,slitD,width,thickness,z_slit2)
        pd = parts.Photodiode(d_pd,z_pd)
        
        
        sub_list_zb = [] #最深到達点保存用
        sub_list_x = []  #PDで計測られた資料表面の光子位置
        sub_list_y = []
        
        try:
            chunk = 50000    #chunk行づつファイルの読み込みを行う（ファイルが大きすぎるため一気に行うと危険）
            csv_obj = pa.read_csv(fname, chunksize = chunk) #CSVファイルの読み込み
            a = np.zeros((chunk,7))                         #読み込んだ後の保存用リスト
            #ファイル読み込み開始
            for r in csv_obj:
                data = list(r.values.flatten())
                for j in range(int(len(data)/7)):
                    if data[7*j] == 0: #エラー箇所の読み飛ばし
                        continue
                    for k in range(7):
                        a[j][k] = data[j*7 + k]
                #フレームデータとして読み込んだdataをlistに直す
                data = list(map(list, a))
                data = [[float(w) for w in v] for v in data]
                
                for i in range (len(a)):
                    if data[i][0] == 0:
                        continue
                    x = data[i][0]; y = data[i][1]; z = 0              #試料表面での光子の位置
                    ux = data[i][2]; uy = data[i][3]; uz = data[i][4]  #試料表面でのベクトル
                    w = data[i][5]; z_b = data[i][6]                   #資料表面での光量、光子最深部到達点
                    
                    if(x!=x)or(y!=y)or(ux!=ux)or(uy!=uy)or(uz!=uz)or(w!=w)or(z_b!=z_b):#エラー（虚数）の除去
                        continue
                    
                    if uz == 0:
                        continue
                    out_count = out_count + 1
                    
                    #スリット１の計算
                    x,y,z = slit_1.hittingPotision(x,y,z,ux,uy,uz,slit_1.front_z)
                    dist = math.sqrt(x**2 + y**2)
                    if (dist <= slit_1.d_in) or (dist >= slit_1.d_out):
                        continue
                    
                    x,y,z = slit_1.hittingPotision(x,y,z,ux,uy,uz,slit_1.back_z)
                    dist = math.sqrt(x**2 + y**2)
                    if (dist <= slit_1.d_in) or (dist >= slit_1.d_out):
                        continue
        
                    #レンズ１での計算
                    x,y,z,ux,uy,uz,leReflect = lens_1.opticalAnalysis(x,y,z,ux,uy,uz)
                    if leReflect:
                        continue
                    
                    #スリット２の計算
                    x,y,z = slit_2.hittingPotision(x,y,z,ux,uy,uz,slit_2.front_z)
                    dist = math.sqrt(x**2 + y**2)
                    if (dist <= slit_2.d_in) or (dist >= slit_2.d_out):
                        continue
    
                    x,y,z = slit_2.hittingPotision(x,y,z,ux,uy,uz,slit_2.back_z)
                    dist = math.sqrt(x**2 + y**2)
                    if (dist <= slit_2.d_in) or (dist >= slit_2.d_out):
                        continue
    
                    #レンズ２の計算
                    x,y,z,ux,uy,uz,leReflect = lens_2.opticalAnalysis(x,y,z,ux,uy,uz)
                    if leReflect:
                        continue
    
                    #フォトダイオードの計算
                    x,y,z = slit_2.hittingPotision(x,y,z,ux,uy,uz,pd.position)
                    dist = math.sqrt(x**2 + y**2)
    
                    if (dist >= pd.r) or (y > 2) or (y < -1):
                        continue
                    #PDで計測された光子の保存
                    pd.count = pd.count + 1
                    pd.record_w = pd.record_w + w
                    sub_list_zb = sub_list_zb + [z_b]
                    #sub_list_x = sub_list_x + [data[i][0]]
                    #sub_list_y = sub_list_y + [data[i][1]]
    
        finally:
            pass
        #file_object.close()
        
        #まとめ
        list_absorbance[count] = pd.absorbance(nPh)
        list_count[count] =  pd.count
        list_w[count] = pd.record_w
        list_zb = list_zb + [sub_list_zb]
        list_x = list_x + [sub_list_x]
        list_y = list_y + [sub_list_y]
        
        rap_time = time.time() - timer_start
        print (rap_time)
        print(count)
        count = count + 1
    
    
    elapsed_time = time.time() - timer_start
    print(elapsed_time)
    print("END")


def showabs(list_distance,list_absorbance):
    plt.plot(list_distance,list_absorbance,".")
    
def showsurf(list_x,list_y):
    plt.plot(list_x[2],list_y[2],".")

def showdepth(step,list_zb,list_distance,):

    #深さとZの関係　標準偏差（1SD）でエラーバー出力
    list_std = np.zeros(step+1)
    list_ave = np.zeros(step+1)
    list_zb = np.array(list_zb)
    for i in range(step+1):
        list_std[i] = np.std(list_zb[i],ddof=1)
        list_ave[i] = sum(list_zb[i])/len(list_zb[i])
    print(list_std)
    plt.plot(list_distance,list_ave,".",color = "r")
    
    plt.errorbar(list_distance,list_ave, list_std,color = "k")
    plt.show()

def saveabs(fname,list_absorbance,list_distance,list_count,list_w):
    #吸光度分布保存用
    fname2 = '160719Ms300Skin0.2nPh10^7TrueOutPut.csv'
    Ms = "300 /cm"
    d_skin = "2 mm"
    coment = "変更なし"
    f = open(fname2,'w')
    
    try:
        
        #f.write("光子数"); f.write(','); f.write(str(nPh)); f.write('\n')
        f.write("μs"); f.write(','); f.write(Ms); f.write('\n')
        f.write("皮膚厚"); f.write(','); f.write(d_skin); f.write('\n')
        f.write("ファイル名"); f.write(','); f.write(fname2); f.write('\n')
        f.write("モンテファイル名"); f.write(','); f.write(fname); f.write('\n')
        f.write("コメント"); f.write(','); f.write(coment); f.write('\n'); f.write('\n')
        f.write("Z"); f.write(','); f.write("absorbance"); f.write(',')
        f.write("count"); f.write(','); f.write("w"); f.write('\n')
        for j in range(len(list_absorbance)):
            dis = list_distance[j]
            abso = round(list_absorbance[j],6)
            con = round(list_count[j],6)
            w_list = round(list_w[j],6)
            f.write(str(dis))
            f.write(',')
            f.write(str(abso))
            f.write(',')
            f.write(str(con))
            f.write(',')
            f.write(str(w_list))
            f.write('\n')
        time.sleep(1.0)
    finally:
        f.close
    print("End")
    
def histdepth(list_zb):

    #ヒストグラム
    df1 = pa.DataFrame(list_zb[0],columns = ['Z = 0'])
    df2 =pa.DataFrame(list_zb[1],columns = ['Z = 1.5'])
    df3 = pa.DataFrame(list_zb[2],columns = ['Z = 3.0'])
    
    df1["Z = 0"].hist(alpha=0.5)
    df2["Z = 1.5"].hist(alpha=0.5)
    df3["Z = 3.0"].hist(alpha=0.5)
    
    plt.xlabel("cm")
    plt.ylabel("Frequency")
    plt.legend(['Z = 0cm','Z = 1.5cm','Z = 3.0cm'])
    plt.show

    """
    #ヒストグラム２
    df1 = pa.DataFrame(list_zb[0],columns = ['Z = 0'])
    df2 =pa.DataFrame(list_zb[1],columns = ['Z = 1.5'])
    df3 = pa.DataFrame(list_zb[2],columns = ['Z = 3.0'])
    #df4 = pa.DataFrame(list_zb[6],columns = ['Z = 9.0'])
    
    df1.plot(alpha = 0.5, color = "blue",kind = 'hist')
    df2.plot(alpha = 0.5, color = "green",kind = 'hist')
    df3.plot(alpha = 0.5, color = "red",kind = 'hist')
    #df4.plot(alpha = 0.5, color = "y",kind = 'hist')
    """

