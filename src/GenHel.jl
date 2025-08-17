# DESCRIPTION
# Open-loop aircraft dynamics.
#
# INPUT
# - dx: state derivative
# - x: current state vector
# - p: current control vector
# - t: time
#
# OUTPUT
# - xdot: system dynamics
#
# ------------------------------------------------------------------------------

function GenHel!(dx,x,p,t)

    # ------------------------------------------------------------------------------
    #
    # ## 1. 初始化与状态划分
    #
    # ------------------------------------------------------------------------------

    # 控制输入混合
    controls = control_mixing(p)

    # 划分状态向量
    # 机身状态 (fuselage states)
    xf=x[1:NFSTATES]
    # 主旋翼状态 (rotor states)
    xr=x[NFSTATES+1:NFSTATES+NRSTATES]
    # 尾桨状态 (tail rotor states)
    xtr=x[NFSTATES+NRSTATES+1:NFSTATES+NRSTATES+NTRSTATES]

    # ------------------------------------------------------------------------------
    #
    # ## 2. 旋翼干扰计算
    #
    # ------------------------------------------------------------------------------

    # 从机体坐标系到桨轴坐标系的转换矩阵
    Tshaft=[cos(ITHSH*D2R)                 0.0            -sin(ITHSH*D2R)
            sin(ITHSH*D2R)*sin(IPHSH*D2R)  cos(IPHSH*D2R)  cos(ITHSH*D2R)*sin(IPHSH*D2R)
            sin(ITHSH*D2R)*cos(IPHSH*D2R) -sin(IPHSH*D2R)  cos(ITHSH*D2R)*cos(IPHSH*D2R)]

    # 提取机体速度
    u=x[1]
    v=x[2]
    w=x[3]

    # 提取纵向挥舞角
    beta1c=x[15]
    A1F=-beta1c*R2D # 转换符号并转为度

    # 提取平均谐波入流
    lambda0=x[29]

    # 计算桨轴坐标系下的速度
    cg=0.0
    Gif=1.0-cg
    temp=Tshaft*[u;v;w]
    us=temp[1]
    vs=temp[2]
    ws=temp[3]
    muxs=us/(OMEGA*R)
    muzs=ws/(OMEGA*R)

    # 计算等效攻角和侧滑角
    chi=R2D*atan2(muxs,abs(lambda0-muzs))+A1F
    alpha=atan2(w,abs(u))*R2D
    beta=atan2(v,sqrt(u^2+w^2))*R2D
    psif=-beta

    # 限制用于查表的值范围
    chitab=min(max(chi,0.0),100.0)
    A1Ftab=min(max(A1F,-6.0),6.0)

    # 通过二维插值查找干扰流因子
    itp1=interpolate((chi_tab, a1f_tab), ekxf_tab', Gridded(Linear()))
    ekxf=itp1(chitab,A1Ftab)
    itp2=interpolate((chi_tab, a1f_tab), ekzf_tab', Gridded(Linear()))
    ekzf=itp2(chitab,A1Ftab)
    itp3=interpolate((chi_tab, a1f_tab), ekxt_tab', Gridded(Linear()))
    ekxt=itp3(chitab,A1Ftab)
    itp4=interpolate((chi_tab, a1f_tab), ekzt_tab', Gridded(Linear()))
    ekzt=itp4(chitab,A1Ftab)

    # 计算旋翼干扰速度和动压损失因子
    rotor_if=zeros(8)
    rotor_if[1]=-Gif*ekzf*lambda0*OMEGA*R  # 干扰流 z 分量
    rotor_if[2]=Gif*ekxf*lambda0*OMEGA*R  # 干扰流 x 分量
    rotor_if[3]=-Gif*ekzt*lambda0*OMEGA*R  # 尾部干扰流 z 分量
    rotor_if[4]=Gif*ekxt*lambda0*OMEGA*R  # 尾部干扰流 x 分量
    rotor_if[5]=table_lookup(alpha_tab,qlossht_tab,alpha) # 平尾动压损失
    rotor_if[6]=table_lookup(psivt_tab,qlossvt_tab,psif)  # 垂尾动压损失
    rotor_if[7]=table_lookup(alphaeps_tab,eps_tab,alpha) # 下洗角
    rotor_if[8]=table_lookup(psisig_tab,sig_tab,psif)   # 侧洗角

    # ------------------------------------------------------------------------------
    #
    # ## 3. 机身气动力计算
    #
    # ------------------------------------------------------------------------------

    # 在配平与线性化过程中，外部速度设为零
    VextF=zeros(3)
    VextTR=zeros(3)
    VextHT=zeros(3)
    VextVT=zeros(3)

    # 机身几何参数 [ft]
    FWT = (FSCG-FSfus)/12
    WWT = (WLCG-WLfus)/12
    BWT = (BLCG-BLfus)/12

    # 阵风速度（当前设为零）
    VXGWF = 0.0
    VYGWF = 0.0
    VZGWF = 0.0

    # NED (北-东-地) 坐标系到机体坐标系的转换
    sphi=sin(xf[7])
    cphi=cos(xf[7])
    sthe=sin(xf[8])
    cthe=cos(xf[8])
    spsi=sin(xf[9])
    cpsi=cos(xf[9])
    TNED2body=[cthe*cpsi                  cthe*spsi                  -sthe
               (sphi*sthe*cpsi-cphi*spsi) (sphi*sthe*spsi+cphi*cpsi)  sphi*cthe
               (cphi*sthe*cpsi+sphi*spsi) (cphi*sthe*spsi-sphi*cpsi)  cphi*cthe]

    # 将外部速度转换到机体坐标系
    Vextb=TNED2body*VextF

    # 考虑旋翼干扰后的机身速度
    VXIWF = rotor_if[2]+Vextb[1]
    VYIWF = 0.0+Vextb[2];
    VZIWF = rotor_if[1]+Vextb[3]

    # 机体坐标系下的速度
    VXB=xf[1]
    VYB=xf[2]
    VZB=xf[3]

    # 计算作用在机身上的总气流速度
    VXWF = VXB+VXGWF+VXIWF
    VYWF = VYB+VYGWF+VYIWF
    VZWF = VZB+VZGWF+VZIWF

    # 计算迎角 (angle of attack) 和侧滑角 (sideslip)
    AVXWF=abs(VXWF);
    if AVXWF < 0.000001
        AVXWF = 0.000001
    end
    ALFWFR = atan(VZWF/AVXWF)
    ALFWF = R2D*ALFWFR
    FVTERM = VXWF^2+VZWF^2
    if FVTERM < 0.000001
        FVTERM = 0.000001
    end
    RVTERM = 1/sqrt(FVTERM)
    CALFWF = VXWF*RVTERM
    SALFWF = VZWF*RVTERM
    BETWFR = atan(VYWF*RVTERM)
    BETWF  = R2D*BETWFR
    SBETWF = sin(BETWFR)
    CBETWF = cos(BETWFR)
    PSIWF = -BETWF

    # 计算动压
    QWF = 0.5*RHO*(VYWF^2 + FVTERM)

    # 通过查表获取机身气动力系数
    AALFWF = abs(ALFWF)
    APSIWF = abs(PSIWF)
    SGNPSI = sign(PSIWF)
    # 迎角相关的系数
    CDA = table_lookup(FUSEAOA,FUSEDA,ALFWF) # 阻力
    CLA = table_lookup(FUSEAOA,FUSELA,ALFWF) # 升力
    CMA = table_lookup(FUSEAOA,FUSEMA,ALFWF) # 俯仰力矩
    # 侧滑角相关的系数
    CLB = table_lookup(FUSEBETA,FUSELB,PSIWF) # 滚转力矩
    CMB = table_lookup(FUSEABETA,FUSEMB,APSIWF)
    CDB = table_lookup(FUSEABETA,FUSEDB,APSIWF)
    CYB = table_lookup(FUSEBETA,FUSEYB,PSIWF) # 侧力
    CRB = table_lookup(FUSEBETA,FUSERB,PSIWF)
    CNB = table_lookup(FUSEBETA,FUSENB,PSIWF) # 偏航力矩

    # 合成总的气动力系数
    CDTOT = CDA + CDB
    CLTOT = CLA + CLB
    CYTOT = CYB
    CRTOT = CRB
    CMTOT = CMA + CMB
    CNTOT = CNB

    # 将气动力从风轴系转换到机体轴系
    Tbw=[CALFWF*CBETWF  CALFWF*SBETWF -SALFWF
         SBETWF        -CBETWF         0.0
         SALFWF*CBETWF  SALFWF*SBETWF  CALFWF]
    Ffw=-QWF*[CDTOT;CYTOT;CLTOT]
    Ff=Tbw*Ffw
    Mf=Tbw*QWF*[CRTOT;-CMTOT;CNTOT]+[0.0 -WWT  BWT; WWT 0.0 -FWT; -BWT FWT 0.0]*Ff
    pfusNED=xf[10:12]+TNED2body'*[FWT;BWT;WWT]

    # ------------------------------------------------------------------------------
    #
    # ## 4. 尾桨气动力计算
    #
    # ------------------------------------------------------------------------------

    # 尾桨模型参数初始化
    VKT = sqrt(u^2+w^2)*fps2kts # 绝对速度 [kts]
    TWSTTR = twistTR/D2R # 尾桨扭转角 [rad]
    OMGRAT = 1 # 主/尾桨转速比
    THETTRC = controls[4] # 尾桨总距

    # 贝利(Bailey)方法所需的大量常数
    C12TR = 0.5
    C13TR = 0.333333333
    C14TR = 0.25
    C15TR = 0.2
    C16TR = 0.166666666
    C23TR = 0.666666666
    C25TR = 0.4
    C43TR = 1.333333333
    C54TR = 1.25
    C58TR = 0.625
    C83TR = 2.666666666
    C89TR = 0.888888888
    XKINF = 4.0/(3.0*pi) # 动态入流常数

    # 机体角速率
    PB=xf[4]
    QB=xf[5]
    RB=xf[6]

    Vextb=TNED2body*VextTR # 外部速度

    # 尾桨几何位置 [ft]
    XTR = -(FSTR-FSCGB)/12.0
    ZTR = -(WLTR-WLCGB)/12.0
    YTR = (BLTR-BLCGB)/12.0

    # 从机体到尾桨的坐标转换（考虑安装角）
    cgam=cos((90-CantTR)*D2R)
    sgam=sin((90-CantTR)*D2R)
    Tcant=[1.0  0.0  0.0
           0.0  cgam sgam
           0.0 -sgam cgam]

    # 更多贝利方法所需的常数
    GTR = a0TR*NB*CHRDTR/(2.0*pi*RTR)
    S0 = a0TR*CHRDTR*RTR^4/IbTR
    S1 = BTLTR/2.0
    S2 = S1*BTLTR
    S3 = S2^2
    S4 = S2/2
    S5 = (BTLTR^3)/3.0
    S6 = DELTTR*TD3TR;
    CBLK = 1.0/(VBVTTR^2)
    B2 = BTLTR^2
    B3 = BTLTR*B2
    B4 = B2^2
    B5 = B2*B3
    B6 = B3^2
    B7 = B3*B4
    B8 = B4^2
    B9 = B4*B5
    B10 = B5^2
    B11 = B5*B6
    B12 = B6^2
    S10 = B2/1296.0
    S11 = C83TR*BTLTR
    S12 = B9/864.0
    S13 = 2.0*B2
    S14 = B10/1080.0
    S15 = C89TR*B2
    S16 = B10/2304.0
    S17 = C43TR*B3
    S18 = B11/1440.0
    S19 = B4/2.0
    S20 = B12/3600.0
    S21 = -B5/108.0
    S22 = -B6/144.0
    S23 = -B7/180.0
    S24 = B2/36.0
    S25 = -C14TR+(1.0/B2)+(C12TR/B4)
    S26 = (B4/162.0)-(B5/81.0)+(B6/144.0)
    S27 = C43TR/BTLTR*(1.0+(1.0/B2))
    S28 = (B5/108.0)-(B6/54.0)+(B7/96.0)
    S29 = 1.0+(1.0/B2)
    S30 = (B6/135.0)*(1.0-(2.0*BTLTR))+(B8/120.0)
    S31 = C14TR+(C89TR/B2)
    S32 = (B6/288.0)-(B7/144.0)+(B8/256.0)
    S33 = C13TR+(C43TR/BTLTR)
    S34 = (B7/180.0)-(B8/90.0)+(B9/160.0)
    S35 = (B8/450.0)-(B9/225.0)+(B10/400.0)
    BLDSTR = NB
    SIGTR = BLDSTR*CHRDTR/(pi*RTR)
    THT1TR = TWSTTR*D2R
    S7 = pi*((OmegaTR*OMGRAT*RTR^2)^2)
    S8 = 1.0/(OmegaTR*RTR*OMGRAT)
    S9 = S7*RTR

    # 计算尾桨处的局部气流速度
    idxtail = 1
    xtail=-(FSTR[idxtail]-FSCG)/12.0
    ytail=(BLTR[idxtail]-BLCG)/12.0
    ztail=-(WLTR[idxtail]-WLCG)/12.0
    qlossfac=sqrt(rotor_if[6])
    ublocal=xf[1]*qlossfac+Vextb[1]+rotor_if[2]+(QB)*ztail-(RB)*ytail
    vblocal=(xf[2]-D2R*xf[1]*rotor_if[8])*qlossfac+Vextb[2]-(PB)*ztail+(RB)*xtail;
    wblocal=(xf[3]-D2R*xf[1]*rotor_if[7])*qlossfac+Vextb[3]+rotor_if[1]+(PB)*ytail-(QB)*xtail

    # 将局部速度转换到尾桨坐标系
    temp=Tcant*[ublocal;vblocal;wblocal]
    ut=temp[1]
    vt=temp[2]
    wt=temp[3]

    # 计算无量纲前进比
    XMUXTR = ut*S8
    XMUYTR = vt*S8
    XMUZTR = wt*S8
    XMU2 = ((XMUXTR^2) + (XMUYTR^2))

    # 计算贝利系数
    BT31 = S2+(0.25*XMU2)
    BT32 = S5+(S1*XMU2)
    BT33 = S3+(S4*XMU2)
    # 用于扭矩计算的额外贝利系数
    XGAM = RHO*S0
    XGAM2 = XGAM^2
    BT41  = S2+(C54TR+(S10*XGAM2))*XMU2
    BT42  = S5+(S11+(S12*XGAM2))*XMU2
    BT43  = S3+(S13+(S14*XGAM2))*XMU2
    BT44  = (S15+(S16*XGAM2))*XMU2
    BT45  = (S17+(S18*XGAM2))*XMU2
    BT46  = (S19+(S20*XGAM2))*XMU2
    BT47  = S21*XMU2
    BT48  = S22*XMU2
    BT49  = S23*XMU2
    BT410 = S24*XMU2
    BT51  = C14TR*(1.0+XMU2)
    BT52  = C13TR
    BT53  = BT51
    BT54  = C15TR+(C16TR*XMU2)
    BT55  = C12TR+((S25+(S26*XGAM2))*XMU2)
    BT56  = C23TR+((S27+(S28*XGAM2))*XMU2)
    BT57  = C12TR+((S29+(S30*XGAM2))*XMU2)
    BT58  = C14TR+((S31+(S32*XGAM2))*XMU2)
    BT59  = C25TR+((S33+(S34*XGAM2))*XMU2)
    BT510 = C16TR+((C58TR+(S35*XGAM2))*XMU2)

    # 垂尾遮挡因子
    BLKTR = BVTTR1
    if VKT < VBVTTR
            BLKTR = ((1.0-BVTTR)*CBLK*(VKT^2))+BVTTR
    end

    # 迭代求解稳态诱导速度和拉力
    TTR = 200
    DWTRSS = xtr[1]
    XLAMTR=XMUZTR-xtr[1]
    THETTR = D2R*(THETTRC)
    iter=0
    err=10
    while abs(err)>1e-6 && iter<200
          # 尾桨桨叶角 [rad]
          THETTR = D2R*(THETTRC-(TTR*S6)+BIASTR)
          # 尾桨下洗和拉力
          DWTRSS = GTR*((XMUZTR*BT31)+(THETTR*BT32)+(D2R*TWSTTR*BT33))/
                   (2.0*sqrt(XMU2+((XMUZTR-DWTRSS)^2))+(GTR*BT31))
          XLAMTR = XMUZTR-DWTRSS
          TTR_new = 2.0*DWTRSS*sqrt(XMU2+(XLAMTR^2))*RHO*S7*BLKTR
          err = TTR_new-TTR
          # 该拉力为稳态版本，之后会被当前入流值覆盖
          TTR = TTR+0.5*err
          iter=iter+1
    end
    if iter>=199
        @printf("warning: TR not converged") # 警告：尾桨未收敛
    end

    # 计算尾桨动态入流
    VT = sqrt(XMU2+(XLAMTR^2))
    CDWTR = XKINF/(OmegaTR*OMGRAT*VT)
    XLAMTR = XMUZTR-xtr[1] # 入流比
    xtrdot=(DWTRSS-xtr[1])/(CDWTR) # 尾桨状态导数

    # 计算尾桨拉力
    CTHTR = GTR*((XLAMTR*BT31)+(THETTR*BT32)+(D2R*TWSTTR*BT33))
    TTR = CTHTR*RHO*S7*BLKTR

    # 计算尾桨扭矩
    ATR = a0TR
    CQTR = (((D2TR*BT55)-(ATR*BT41))*(XLAMTR^2)+((D2TR*BT56)-
           (ATR*BT42))*THETTR*XLAMTR+((D2TR*BT57)-(ATR*BT43))*THT1TR*XLAMTR+
           ((D2TR*BT58)-(ATR*BT44))*(THETTR^2)+((D2TR*BT59)-
           (ATR*BT45))*THT1TR*THETTR+((D2TR*BT510)-(ATR*BT46))*(THT1TR^2)+
           (D0TR*BT51 )+D1TR*((BT52*XLAMTR)+(BT53*THETTR)+
           (BT54*THT1TR)))*0.5*SIGTR
    Qtr = CQTR*RHO*S9

    # 计算尾桨阻力
    DTR=0.5*CDTR*RHO*(ut^2)
    # 计算尾桨功率 [hp]
    HPTR = Qtr*OmegaTR*OMGRAT*FtLb_s2Hp

    # 将尾桨的力和力矩转换到机体坐标系
    Ftr=Tcant'*[-DTR;0.0;-TTR]
    Mtr_local=Tcant'*[0.0;0.0;-Qtr]
    Mtr=[0.0  -ZTR  YTR
         ZTR   0.0 -XTR
         -YTR  XTR  0.0]*Ftr+Mtr_local
    ptrNED = xf[10:12]+TNED2body'*[XTR;YTR;ZTR]

    # ------------------------------------------------------------------------------
    #
    # ## 5. 平尾 (Horizontal Tail) 气动力计算
    #
    # ------------------------------------------------------------------------------

    idxtail=1 # 平尾索引
    # 根据飞行速度和俯仰速率计算安定面偏角
    Veq=sqrt( (x[1])^2 + (x[2])^2 + (x[3])^2 )/1.688*sqrt(RHO/rhoSLSTD)
    ayg=(dx[1]-G*sin(x[7])*cos(x[8])+x[6]*x[1]-x[4]*x[3])/G
    StabInc=stabSched(p[3],Veq,ayg,x[5],STABSET)

    # 平尾几何位置
    xtail=-(FSTAIL[idxtail]-FSCG)/12.0
    ytail=(BLTAIL[idxtail]-BLCG)/12.0
    ztail=-(WLTAIL[idxtail]-WLCG)/12.0
    Vextb=TNED2body*VextHT

    # 计算平尾处的局部气流速度（考虑旋翼下洗和动压损失）
    qlossfac=sqrt(rotor_if[5])
    epsilon=rotor_if[7]
    sigma=0.0
    uHT=(xf[1])*qlossfac+Vextb[1]+rotor_if[4]
    vHT=(xf[2]-D2R*xf[1]*sigma)+Vextb[2]
    wHT=(xf[3]-D2R*xf[1]*epsilon)*qlossfac+Vextb[3]+rotor_if[3]
    ublocal=uHT+(QB)*ztail-(RB)*ytail
    vblocal=vHT-(PB)*ztail+(RB)*xtail
    wblocal=wHT+(PB)*ytail-(QB)*xtail

    # 速度坐标转换
    if (CLELEV[idxtail] == 0.0)
        cthe=cos((ITAIL[idxtail]+StabInc)*D2R)
        sthe=sin((ITAIL[idxtail]+StabInc)*D2R)
    else
        cthe=cos(ITAIL[idxtail]*D2R)
        sthe=sin(ITAIL[idxtail]*D2R)
    end
    cphi=cos(PHITAIL[idxtail]*D2R)
    sphi=sin(PHITAIL[idxtail]*D2R)
    T=[cthe  sthe*sphi -sthe*cphi
       0     cphi       sphi
       sthe -sphi*cthe  cphi*cthe]
    temp=T*[ublocal;vblocal;wblocal]
    ut=temp[1]
    vt=temp[2]
    wt=temp[3]
    Vtot2=ut^2+vt^2+wt^2
    Vtot=sqrt(Vtot2)
    Vxz=sqrt(ut^2+wt^2)
    qt=0.5*RHO*Vtot2

    # 计算平尾的迎角和侧滑角
    alphat=atan2(wt,ut)*R2D
    sa=wt/Vxz
    ca=ut/Vxz
    sb=vt/Vtot
    cb=Vxz/Vtot

    # 查表计算升力和阻力系数
    CLt=table_lookup(ALTAIL[idxtail,:],CLTAIL[idxtail,:],alphat)+CLELEV[idxtail]*StabInc
    CDt=table_lookup(ALTAIL[idxtail,:],CDTAIL[idxtail,:],alphat)
    Lift=qt*STAIL[idxtail]*CLt
    Drag=qt*STAIL[idxtail]*CDt

    # 将力和力矩转换回机体坐标系
    Tw2t=[cb*ca -sb  -cb*sa
          sb*ca  cb  -sb*sa
          sa     0.0  ca]
    temp=T'*Tw2t*[-Drag;0.0;-Lift]
    Fht=temp
    Mht=[ 0.0   -ztail  ytail
         ztail  0.0   -xtail
        -ytail  xtail  0.0]*Fht
    phtNED=xf[10:12]+TNED2body'*[xtail;ytail;ztail]

    # ------------------------------------------------------------------------------
    #
    # ## 6. 垂尾 (Vertical Tail) 气动力计算
    #
    # ------------------------------------------------------------------------------

    idxtail=2 # 垂尾索引
    # 垂尾几何位置
    xtail=-(FSTAIL[idxtail]-FSCG)/12.0
    ytail=(BLTAIL[idxtail]-BLCG)/12.0
    ztail=-(WLTAIL[idxtail]-WLCG)/12.0
    Vextb=TNED2body*VextVT

    # 计算垂尾处的局部气流速度（考虑旋翼侧洗和动压损失）
    qlossfac=sqrt(rotor_if[6])
    epsilon=0.0
    sigma=rotor_if[8]
    uHT=(xf[1])*qlossfac+Vextb[1]+rotor_if[4]
    vHT=(xf[2]-D2R*xf[1]*sigma)+Vextb[2]
    wHT=(xf[3]-D2R*xf[1]*epsilon)*qlossfac+Vextb[3]+rotor_if[3]
    ublocal=uHT+(QB)*ztail-(RB)*ytail
    vblocal=vHT-(PB)*ztail+(RB)*xtail
    wblocal=wHT+(PB)*ytail-(QB)*xtail

    # 速度坐标转换
    if (CLELEV[idxtail] == 0.0)
        # 垂尾安装角设为0
        cthe=cos((ITAIL[idxtail]+0.0)*D2R)
        sthe=sin((ITAIL[idxtail]+0.0)*D2R)
    else
        cthe=cos(ITAIL[idxtail]*D2R)
        sthe=sin(ITAIL[idxtail]*D2R)
    end
    cphi=cos(PHITAIL[idxtail]*D2R)
    sphi=sin(PHITAIL[idxtail]*D2R)
    T=[cthe  sthe*sphi -sthe*cphi
       0     cphi       sphi
       sthe -sphi*cthe  cphi*cthe]
    temp=T*[ublocal;vblocal;wblocal]
    ut=temp[1]
    vt=temp[2]
    wt=temp[3]
    Vtot2=ut^2+vt^2+wt^2
    Vtot=sqrt(Vtot2)
    Vxz=sqrt(ut^2+wt^2)
    qt=0.5*RHO*Vtot2

    # 计算垂尾的迎角和侧滑角
    alphat=atan2(wt,ut)*R2D
    sa=wt/Vxz
    ca=ut/Vxz
    sb=vt/Vtot
    cb=Vxz/Vtot

    # 查表计算升力和阻力系数
    CLt=table_lookup(ALTAIL[idxtail,:],CLTAIL[idxtail,:],alphat)+CLELEV[idxtail]*0.0
    CDt=table_lookup(ALTAIL[idxtail,:],CDTAIL[idxtail,:],alphat)
    Lift=qt*STAIL[idxtail]*CLt
    Drag=qt*STAIL[idxtail]*CDt

    # 将力和力矩转换回机体坐标系
    Tw2t=[cb*ca -sb  -cb*sa
          sb*ca  cb  -sb*sa
          sa     0.0  ca]
    temp=T'*Tw2t*[-Drag;0.0;-Lift]
    Fvt=temp
    Mvt=[ 0.0   -ztail  ytail
         ztail  0.0   -xtail
        -ytail  xtail  0.0]*Fvt
    pvtNED=xf[10:12]+TNED2body'*[xtail;ytail;ztail]

    # ------------------------------------------------------------------------------
    #
    # ## 7. 主旋翼 (Main Rotor) 气动力与动力学计算
    #
    # ------------------------------------------------------------------------------

    # --- 7.1 输入与初始化 ---

    # 定义传递给主旋翼模块的输入向量 ur
    ur=zeros(23+NB*NSEG*3+1)
    ur[1]=controls[1]*pi/180 # 横向周期变距
    ur[2]=controls[2]*pi/180 # 纵向周期变距
    ur[3]=controls[3]*pi/180 # 总距
    ur[4:15]=xf[1:12]         # 机身状态
    ur[16:21]=dx[1:6]          # 机身加速度
    ur[22]=OMEGA               # 旋翼名义转速
    ur[23]=0                   # 旋翼角加速度 (当前为0)
    ur[24:23+NB*NSEG*3+1]=zeros(NB*NSEG*3+1,1) # 外部流场速度 (阵风)

    if ur[23]!=0.0
        soptval=1
    end

    # 根据桨根和桨尖的扭转角，插值计算每个桨段的扭转角
    itp = interpolate((TWISTTABR,),TWISTTABTHET, Gridded(Linear()))
    twist = itp(RSEG)

    # --- 7.2 桨叶运动学 ---

    # 从状态向量中提取主旋翼状态
    beta1c=xr[3]  # 纵向挥舞角
    lambda0=xr[17] # 平均诱导流
    lambda1s=xr[18] # 侧向谐波诱导流
    lambda1c=xr[19] # 纵向谐波诱导流
    psi1=xr[20]   # 旋翼方位角

    # 计算每个桨叶的方位角
    psi=psi1*ones(4)+collect(0:(NB-1))*(2*pi/NB)
    for i=1:NB
        psi[i]=mod(psi[i]',2*pi)
    end

    # 计算方位角的三角函数
    cpsi = zeros(4)
    spsi = zeros(4)
    for i=1:NB
        cpsi[i]=cos(psi[i])
        spsi[i]=sin(psi[i])
    end

    # 桨叶运动转换：从多桨叶坐标系 (MBC) 到独立桨叶坐标系 (IBC)
    LMBC2IBC=[ones(4) (-1).^collect(1:NB) cpsi spsi]
    LMBC2IBCdot=OMEGA*[zeros(4,2) -spsi cpsi]
    beta=LMBC2IBC*xr[1:NB] # 挥舞角
    beta_dot=LMBC2IBC*xr[NB+1:2*NB]+LMBC2IBCdot*xr[1:NB] # 挥舞角速度
    zeta=LMBC2IBC*xr[2*NB+1:3*NB] # 摆振角
    zeta_dot=LMBC2IBC*xr[3*NB+1:4*NB]+LMBC2IBCdot*xr[2*NB+1:3*NB] # 摆振角速度

    # 桨叶运动转换：从独立桨叶坐标系 (IBC) 到多桨叶坐标系 (MBC)
    LIBC2MBC=1/NB*[ones(4)'
                   ((-1).^collect(1:NB))'
                   2*cpsi'
                   2*spsi']
    LIBC2MBCdot=OMEGA/NB*[zeros(2,4)
                           -2*spsi'
                           2*cpsi']
    LIBC2MBCddot=(OMEGA*OMEGA)/NB*[zeros(2,4)
                                   -2*cpsi'
                                   -2*spsi']

    # 计算桨叶角度的三角函数
    cbeta = zeros(4)
    sbeta = zeros(4)
    czeta = zeros(4)
    szeta = zeros(4)
    cpsizeta = zeros(4)
    spsizeta = zeros(4)
    for i=1:NB
        cbeta[i]=cos(beta[i])
        sbeta[i]=sin(beta[i])
        czeta[i]=cos(zeta[i])
        szeta[i]=sin(zeta[i])
        cpsizeta[i]=cos(psi[i]+zeta[i])
        spsizeta[i]=sin(psi[i]+zeta[i])
    end
    # --- 7.3 桨叶局部速度计算 ---

    # 从输入向量 ur 中提取各物理量
    theta1c=ur[1] # 横向周期变距
    theta1s=ur[2] # 纵向周期变距
    theta0=ur[3]  # 总距
    Vxb=ur[4]     # 机体x轴速度
    Vyb=ur[5]     # 机体y轴速度
    Vzb=ur[6]     # 机体z轴速度
    p=ur[7]       # 滚转角速率
    q=ur[8]       # 俯仰角速率
    r=ur[9]       # 偏航角速率
    sphi=sin(ur[10])
    cphi=cos(ur[10])
    sthe=sin(ur[11])
    cthe=cos(ur[11])
    spsiE=sin(ur[12])
    cpsiE=cos(ur[12])
    Tbody2NED=[cthe*cpsiE (sphi*sthe*cpsiE-cphi*spsiE) (cphi*sthe*cpsiE+sphi*spsiE)
               cthe*spsiE (sphi*sthe*spsiE+cphi*cpsiE) (cphi*sthe*spsiE-sphi*cpsiE)
               -sthe                  sphi*cthe                 cphi*cthe]
    xcg=ur[13]    # 重心位置 x
    ycg=ur[14]    # 重心位置 y
    zcg=ur[15]    # 重心位置 z
    Vxbdot=ur[16] # 机体x轴加速度
    Vybdot=ur[17] # 机体y轴加速度
    Vzbdot=ur[18] # 机体z轴加速度
    pdot=ur[19]   # 滚转角加速度
    qdot=ur[20]   # 俯仰角加速度
    rdot=ur[21]   # 偏航角加速度
    Omega_dot=ur[23] # 旋翼角加速度
    VgNED=reshape(ur[24:23+NB*NSEG*3],3,NB*NSEG) # 阵风速度
    InflowFade=ur[23+NB*NSEG*3+1] # 入流模型淡出增益

    # 将阵风速度转换到桨轴坐标系
    Vgs=Tshaft*Tbody2NED'*VgNED;
    Vgxs=reshape((Vgs[1,:])',NSEG,NB)'
    Vgys=reshape((Vgs[2,:])',NSEG,NB)'
    Vgzs=reshape((Vgs[3,:])',NSEG,NB)'

    # 计算桨毂处的加速度
    gx=G*sthe
    gy=-G*sphi*cthe
    gz=-G*cphi*cthe
    Xh=(FSCGB-FSMR)/12.0 # 从重心到桨毂的力臂
    Yh=(BLCGB-BLMR)/12.0
    Zh=(WLCGB-WLMR)/12.0
    Vxhdot=Vxbdot-r*Vyb+q*Vzb-Xh*(q^2+r^2)+Yh*(p*q-rdot)+Zh*(p*r+qdot)+gx
    Vyhdot=Vybdot-p*Vzb+r*Vxb+Xh*(p*q+rdot)-Yh*(p^2+r^2)+Zh*(q*r-pdot)+gy
    Vzhdot=Vzbdot+p*Vyb-q*Vxb+Xh*(p*r-qdot)+Yh*(q*r+pdot)-Zh*(p^2+q^2)+gz

    # 计算无量纲前进比
    muxh=(Vxb+q*Zh-r*Yh)/(OMEGA*R)
    muyh=(Vyb+r*Xh-p*Zh)/(OMEGA*R)
    muzh=(Vzb-q*Xh+p*Yh)/(OMEGA*R)
    muxhdot=(Vxbdot+qdot*Zh-rdot*Yh)/(OMEGA*R) # 前进比导数
    muyhdot=(Vybdot+rdot*Xh-pdot*Zh)/(OMEGA*R)
    muzhdot=(Vzbdot-qdot*Xh+pdot*Yh)/(OMEGA*R)

    # 将速度和加速度转换到桨轴坐标系
    temp=Tshaft*[Vxhdot;Vyhdot;Vzhdot]
    Vxsdot=temp[1]
    Vysdot=temp[2]
    Vzsdot=temp[3]
    temp=Tshaft*[pdot;qdot;rdot]
    psdot=temp[1]
    qsdot=temp[2]
    rsdot=temp[3]
    rsdotmomd=rsdot-Omega_dot
    temp=Tshaft*[muxh;muyh;muzh]
    muxs=temp[1]
    muys=temp[2]
    muzs=temp[3]
    mu=sqrt(muxs^2+muys^2)
    betawind = atan2(muys, muxs)
    temp=Tshaft*[muxhdot;muyhdot;muzhdot]
    muxsdot=temp[1]
    muysdot=temp[2]
    muzsdot=temp[3]
    betawind_dot=(muxs*muysdot-muys*muxsdot)/(mu^2)
    temp=Tshaft*[p;q;r]
    ps=temp[1]
    qs=temp[2]
    rs=temp[3]
    rminusO=rs-OMEGA

    # 计算每个桨叶的桨距角及其变化率
    thetab = zeros(4)
    thetab_dot = zeros(4)
    for i=1:NB
        thetab[i]=theta0+theta1c*cos(psi[i]+DELSP)+theta1s*sin(psi[i]+DELSP)-
                  tan(DELTA3)*beta[i]
        thetab_dot[i]=-theta1c*OMEGA*sin(psi[i]+DELSP)+theta1s*OMEGA*
                      cos(psi[i]+DELSP)
    end

    # 计算动态扭转分量
    Veq=sqrt(RHO/rhoSLSTD)*sqrt(Vxb^2+Vyb^2+Vzb^2)/1.688
    KVDT=min(max(1.4e-4-(4.4e-6*Veq),-0.00052),-0.0003)
    theta_DTtip=KVDT*xr[21]
    theta_DT=vec(ones(4,1))*vec(DynTwistMode)'*theta_DTtip*pi/180.0
    # 合成每个桨段的总桨距角
    theta=vec(thetab)*ones(1,NSEG)+ones(NB,1)*vec(twist)'*pi/180+theta_DT

    # 将阵风速度转换到桨叶坐标系
    Utg=1/(OMEGA*R)*( (spsizeta*ones(1,NSEG)).*Vgxs + (cpsizeta*ones(1,NSEG)).*Vgys)
    Urg=1/(OMEGA*R)*( ((-cpsizeta.*cbeta)*ones(1,NSEG)).*Vgxs +
        ((spsizeta.*cbeta)*ones(1,NSEG)).*Vgys - (sbeta*ones(1,NSEG)).*Vgzs )
    Upg=1/(OMEGA*R)*( ((-cpsizeta.*sbeta)*ones(1,NSEG)).*Vgxs +
        ((spsizeta.*sbeta)*ones(1,NSEG)).*Vgys + (cbeta*ones(1,NSEG)).*Vgzs )

    # 计算 Pitt-Peters 入流模型引起的诱导速度
    Upd_pp=(-lambda0*vec(cbeta)*ones(1,NSEG)-lambda1c*(e*vec(cbeta).*vec(cpsi)*ones(1,NSEG)+
           vec(cbeta).*vec(cpsizeta)*vec(RSEG)')-lambda1s*(e*vec(cbeta).*vec(spsi)*ones(1,NSEG)+vec(cbeta).*
           vec(spsizeta)*vec(RSEG)') )*(1-InflowFade)
    Urd_pp=(-lambda0*vec(sbeta)*ones(1,NSEG)-lambda1c*(e*vec(sbeta).*vec(cpsi)*ones(1,NSEG)+
           vec(sbeta).*vec(cpsizeta)*vec(RSEG)')-lambda1s*(e*vec(sbeta).*vec(spsi)*ones(1,NSEG)+vec(sbeta).*
           vec(spsizeta)*vec(RSEG)') )*(1-InflowFade)
    Upd_ph = zeros(NB,NSEG)
    Urd_ph = zeros(NB,NSEG)

    # 计算每个桨叶剖面上的三个方向的局部气流速度 (Up, Ut, Ur)
    # Up: 垂直于桨盘平面的速度 (挥舞方向)
    # Ut: 桨叶剖面切向速度 (旋转方向)
    # Ur: 桨叶径向速度
    Up=(-muxs*vec(sbeta).*vec(cpsizeta)+muys*vec(sbeta).*vec(spsizeta)+muzs*vec(cbeta))*ones(1,NSEG)+
       (e/OMEGA)*(vec(cbeta).*(qs*vec(cpsi)+ps*vec(spsi))-vec(sbeta).*vec(szeta)*rminusO)*ones(1,NSEG)+
       (-vec(beta_dot)+qs*vec(cpsizeta)+ps*vec(spsizeta))*vec(RSEG)'/OMEGA+Upd_pp+Upg+Upd_ph
    Ut=(muxs*vec(spsizeta)+muys*vec(cpsizeta)-(e/OMEGA)*vec(czeta)*rminusO)*ones(1,NSEG)+
       ((vec(zeta_dot)-ones(4)*rminusO).*vec(cbeta)+(ps*vec(cpsizeta)-qs*vec(spsizeta)).*vec(sbeta))*vec(RSEG)'/OMEGA+
       Utg
    Ur=(muxs*vec(cbeta).*vec(cpsizeta)-muys*vec(cbeta).*vec(spsizeta)+muzs*vec(sbeta))*ones(1,NSEG)+
       (e/OMEGA)*(vec(sbeta).*(qs*vec(cpsi)+ps*vec(spsi))+vec(cbeta).*vec(szeta)*rminusO)*ones(1,NSEG)+
       Urd_pp+Urg+Urd_ph
    Utot=(Up.^2+Ut.^2+Ur.^2).^(1/2) # 总速度

    # --- 7.4 剖面气动力系数计算 ---

    # 计算每个桨段的局部迎角、马赫数和偏流角
    cosgam=abs.(Ut)./((Ut.^2+Ur.^2).^(1/2))
    acosgam=abs.(cosgam)
    alphay=(180.0/pi)*atan2.((Ut.*tan.(theta)+Up).*acosgam,(Ut-Up.*tan.(theta).*
           cosgam.*cosgam))
    Mach=((Ut.^2+Up.^2).^(1/2))*OMEGA*R/VSOUND
    Machtab=max.(min.(Mach,1.0),0.0) # 限制马赫数范围

    # GenHel 旋翼模型中的迎角转换，用于处理大迎角和反流区
    atrans0=cosgam.*alphay
    atranslow=(alphay+ones(4,10)*180.0).*cosgam-ones(4,10)*180.0
    atranshigh=ones(4,10)*180.0+(alphay-ones(4,10)*180.0).*cosgam
    atrans1=((1.0-ACL1/90.0)*atrans0-ACL1*(ones(4,10)-cosgam))./(-ACL1/90.0*ones(4,10)+cosgam)
    atrans2=((1.0-ACL2/90.0)*atranshigh+ACL2*(ones(4,10)-cosgam))./(2.0*ones(4,10)-ACL2/90.0*ones(4,10)-cosgam)
    atrans3=((1.0+ACL3/90.0)*atrans0-ACL3*(ones(4,10)-cosgam))./(ACL3/90.0*ones(4,10)+cosgam)
    atrans4=((1.0+ACL4/90.0)*atranslow+ACL4*(ones(4,10)-cosgam))./(2.0*ones(4,10)+ACL4/90.0*ones(4,10)-cosgam)
    atrans=atrans0.*(atrans0.<=ACL1).*(atrans0.>=ACL3).*(abs.(alphay).<90.0)+
           atrans1.*(atrans0.>ACL1).*(alphay.<90.0)+(atrans2.*(atranshigh.<ACL2)+
           atranshigh.*(atranshigh.>=ACL2)).*(alphay.>=90.0)+atrans3.*
           (atrans0.<ACL3).*(alphay.>=-90.0)+(atrans4.*(atranslow.>ACL4)+
           atranslow.*(atranslow.<=ACL4)).*(alphay.<-90.0)
    atransBtab=max.(min.(atrans,32.0),-32.0)
    alphayBtab=max.(min.(alphay,32.0),-32.0)

    # 查表获取升力系数 (CL) 和阻力系数 (CD)
    # 对小迎角范围使用二维插值(迎角和马赫数)，对大迎角范围使用一维插值
    CL=interpp1(vec(AOAUTAB),vec(CLR0UTAB),atrans).*(abs.(atrans).>32.0)+
       interpp2(vec(AOABTAB),vec(MACHTAB),CLR0BTAB,atransBtab,Machtab).*
       (abs.(atrans).<=32.0)
    CD=interpp1(vec(AOAUTAB),vec(CDR0UTAB),alphay).*(abs.(alphay).>32.0)+
       interpp2(vec(AOABTAB),vec(MACHTAB),CDR0BTAB,alphayBtab,Machtab).*
       (abs.(alphay).<=32.0)+ones(4,10)*DCDMR
    # --- 7.5 剖面气动力和力矩计算 ---

    # Peters-He 入流模型专用的升力数组
    BLSEGLIFT=((0.5*(ones(NB,1)*vec(CHORD)'))./R).*CL.*sign.(Ut).*Utot.*
              sqrt.(Ut.^2+Ur.^2)

    # 计算每个桨段上的气动力
    Fp=(0.5*RHO*OMEGA^2*R^3)*(ones(NB,1)*vec(CHORD)').*Utot.* # 挥舞方向力
       (CL.*Ut./acosgam+CD.*Up).*(ones(NB,1)*vec(delseg)')
    Ft=(0.5*RHO*OMEGA^2*R^3)*(ones(NB,1)*vec(CHORD)').*Utot.* # 切向力 (阻力)
       (CD.*Ut-CL.*Up.*acosgam).*(ones(NB,1)*vec(delseg)')
    Fr=(0.5*RHO*OMEGA^2*R^3)*(ones(NB,1)*vec(CHORD)').*Utot.* # 径向力
       (CD-CL.*Up.*acosgam./Ut).*Ur.*(ones(NB,1)*vec(delseg)')

    # 对每个桨叶的剪力进行积分
    Fpb=zeros(4,1)
    Ftb=zeros(4,1)
    Frb=zeros(4,1)
    for i=1:NB
        Fpb[i]=sum(Fp[i,:])
        Ftb[i]=sum(Ft[i,:])
        Frb[i]=sum(Fr[i,:])
    end

    # 计算用于动态扭转模型的总桨叶力
    Fp_DT=sum(sqrt.(Fpb.^2+Ftb.^2))/NB

    # 对每个挥舞铰的航空力矩进行积分
    Mf_aero=zeros(4,1) # 挥舞力矩
    Ml_aero=zeros(4,1) # 摆振力矩
    for i=1:NB
        Mf_aero[i]=R*sum(rseg.*Fp[i,:])
        Ml_aero[i]=R*sum(rseg.*Ft[i,:])
    end

    # 桨毂处的气动力矩
    Lha=-sum(Mf_aero.*spsizeta) # 滚转力矩
    Mha=-sum(Mf_aero.*cpsizeta) # 俯仰力矩

    # 桨毂处的气动剪力
    Fxa=Frb.*cbeta.*szeta-Ftb.*czeta-Fpb.*sbeta.*szeta
    Fya=Frb.*cbeta.*czeta+Ftb.*szeta-Fpb.*sbeta.*czeta
    Fza=-Frb.*sbeta-Fpb.*cbeta

    # 总气动拉力
    Tha=-sum(Fza)

    # 计算用于 Pitt-Peters 模型的气动力和力矩系数
    CTa=Tha/(pi*RHO*OMEGA^2*R^4) # 拉力系数
    CLa=Lha/(pi*RHO*OMEGA^2*R^5) # 滚转力矩系数
    CMa=Mha/(pi*RHO*OMEGA^2*R^5) # 俯仰力矩系数
    # --- 7.6 桨叶动力学与惯性力计算 ---

    # 摆振阻尼器运动学
    thetald=thetab-THETALDGEO*ones(4)
    xld=ALD*sbeta+ones(4)*CLD+(BLD*cos.(zeta+ones(4)*ZETA0)+DLD*
        sin.(zeta+ones(4)*ZETA0)).*cbeta
    yld=-RLD*cos.(thetald)-BLD*sin.(zeta+ones(4)*ZETA0)+DLD*
        cos.(zeta+ones(4)*ZETA0)
    zld=ALD*cbeta-RLD*sin.(thetald)-(BLD*cos.(zeta+ones(4)*ZETA0)+
        DLD*sin.(zeta+ones(4)*ZETA0)).*sbeta
    daxld=sqrt.(xld.^2+yld.^2+zld.^2)
    xld_dot = ALD*beta_dot.*cbeta-beta_dot.*sbeta.*(BLD*cos.(zeta+ones(4)
              *ZETA0)+DLD*sin.(zeta+ones(4)*ZETA0))+(-BLD*sin.(zeta+ones(4)*
              ZETA0)+DLD*cos.(zeta+ones(4)*ZETA0)).*zeta_dot.*cbeta
    yld_dot = RLD*sin.(thetald).*thetab_dot-(BLD*cos.(zeta+ones(4)*ZETA0)+DLD*
              cos.(zeta+ones(4)*ZETA0)).*zeta_dot
    zld_dot = -ALD*sbeta.*beta_dot-RLD*cos.(thetald).*thetab_dot -cbeta.*
              beta_dot.*(BLD*cos.(zeta+ones(4)*ZETA0)+DLD*sin.(zeta+ones(4)*
              ZETA0))-(-BLD*sin.(zeta+ones(4)*ZETA0)+DLD*cos.(zeta+ones(4)*
              ZETA0)).*zeta_dot.*sbeta
    raxld = (xld.*xld_dot+yld.*yld_dot+zld.*zld_dot)./daxld

    # 线性摆振阻尼器产生的力
    Fld=CLAG*raxld

    # 摆振阻尼器产生的挥舞和摆振力矩
    Mf_ld=-Fld.*((zld*CLD)+(xld.*sin.(thetald)*RLD))./(12.0*daxld)
    Ml_ld=-Fld.*(RLD*cos.(thetald).*(xld.*cbeta-zld.*sbeta)+yld.*
          (CLD*cbeta+RLD*sin.(thetald).*sbeta))./(12.0*daxld)-
          KLAG*(zeta+ones(4)*ZETA0)

    # 虚拟的挥舞弹簧/阻尼器力矩
    Mf_fd=-KFLAP*beta-CFLAP*beta_dot

    # 计算挥舞和摆振角加速度 (beta_ddot, zeta_ddot)
    beta_ddot = zeros(NB)
    zeta_ddot = zeros(NB)
    for i=1:NB
        beta_ddot[i]=(MBETA/IBETA)*(cbeta[i].*(Vzsdot+HOFFSET*(2*OMEGA*(ps*cpsi[i]-qs*spsi[i])+
              psdot*spsi[i]+qsdot*cpsi[i]))+sbeta[i].*czeta[i].*(Vysdot*spsi[i]-Vxsdot*cpsi[i]-
              HOFFSET*rminusO^2)+sbeta[i].*szeta[i].*(Vysdot*spsi[i]-Vxsdot*cpsi[i]-
              HOFFSET*rsdotmomd))+(cbeta[i].^2).*(czeta[i].*(psdot*spsi[i]+qsdot*
              cpsi[i]-2*(zeta_dot[i]+OMEGA).*(qs*spsi[i]-ps*cpsi[i]))-szeta[i].*(2*(OMEGA+
              zeta_dot[i]).*(ps*spsi[i]+qs*cpsi[i])+(qsdot*spsi[i]-psdot*cpsi[i])))+
              cbeta[i].*sbeta[i].*(2*zeta_dot[i]*rminusO-rminusO^2-zeta_dot[i].^2)+
              (Mf_aero[i]+Mf_ld[i]+Mf_fd[i])/IBETA
        zeta_ddot[i]=(MBETA./(IBETA*cbeta[i])).*(szeta[i].*(Vysdot*spsi[i]-Vxsdot*cpsi[i]-
                  HOFFSET*rminusO^2)-czeta[i].*(Vxsdot*spsi[i]+Vysdot*cpsi[i]))+rsdotmomd+
                  (sbeta[i]./cbeta[i]).*(2*beta_dot[i].*(OMEGA+zeta_dot[i]-rs)+qsdot*spsizeta[i]-
                  psdot*cpsizeta[i])+2*beta_dot[i].*(czeta[i].*(qs*spsi[i]-ps*cpsi[i])+szeta[i].*
                  (ps*spsi[i]+qs*cpsi[i]))+(Ml_ld[i])./(IBETA*cbeta[i].^2)-
                  Ml_aero[i]./(IBETA*cbeta[i])
    end

    # 计算每个铰链处的惯性剪力
    Fxi = zeros(NB)
    Fyi = zeros(NB)
    Fzi = zeros(NB)
    for i=1:NB
        Fxi[i]=MBETA*(cbeta[i].*czeta[i].*(rsdotmomd-zeta_ddot[i])+2*sbeta[i].*czeta[i].*(zeta_dot[i].*
            beta_dot[i]-rminusO*beta_dot[i])+cbeta[i].*szeta[i].*(zeta_dot[i].^2+beta_dot[i].^2-2*
            rminusO*zeta_dot[i]+rminusO^2)+2*beta_dot[i].*cbeta[i].*(ps*cpsi[i]-qs*spsi[i])+
            sbeta[i].*szeta[i].*beta_ddot[i])-WBLADE/G*(Vxsdot*spsi[i]+Vysdot*cpsi[i])
        Fyi[i]=MBETA*(cbeta[i].*czeta[i].*(zeta_dot[i].^2+beta_dot[i].^2-2*rminusO*zeta_dot[i]+
            rminusO^2)+sbeta[i].*czeta[i].*beta_ddot[i]+cbeta[i].*szeta[i].*zeta_ddot[i]-2*beta_dot[i].*
            cbeta[i].*(ps*spsi[i]+qs*cpsi[i])+WBLADE*HOFFSET/(G*MBETA)*rminusO^2)+WBLADE/G*
            (Vxsdot*cpsi[i]-Vysdot*spsi[i])
        Fzi[i]=MBETA*(beta_ddot[i].*cbeta[i]-(beta_dot[i].^2).*sbeta[i] + 2*sbeta[i].*czeta[i].*
            beta_dot[i].*(ps*spsi[i]+qs*cpsi[i])+cbeta[i].*szeta[i].*(2*(OMEGA+zeta_dot[i]).*
            (ps*spsi[i]+qs*cpsi[i])+qsdot*spsi[i]-psdot*cpsi[i])-cbeta[i].*czeta[i].*(2*(OMEGA+
            zeta_dot[i]).*(ps*cpsi[i]-qs*spsi[i])+psdot*spsi[i]+qsdot*cpsi[i])+WBLADE*HOFFSET/
            (G*MBETA)*((2*OMEGA*(ps*cpsi[i]-qs*spsi[i])+psdot*spsi[i]+qsdot*cpsi[i])))-
            WBLADE/G*Vzsdot
    end

    # 计算总剪力 (惯性力 + 气动力)
    Fxt=Fxi+Fxa
    Fyt=Fyi+Fya
    Fzt=Fzi+Fza

    # 计算桨毂处的总力和力矩
    Thrust=-sum(Fzt) # 拉力
    Hforce=sum(Fyt.*cpsi-Fxt.*spsi) # H-力 (沿y轴)
    Jforce=-sum(Fxt.*cpsi+Fyt.*spsi) # J-力 (沿x轴)
    Lhub=sum(HOFFSET*Fzt.*spsi+Mf_ld.*spsizeta) # 滚转力矩
    Mhub=sum(HOFFSET*Fzt.*cpsi+Mf_ld.*cpsizeta) # 俯仰力矩
    Qhub=-sum(HOFFSET*Fxt-Ml_ld) # 扭矩
    # --- 7.7 Pitt-Peters 动态入流模型 ---

    # 迭代求解准稳态入流 lambda0_qs
    lambda0_qs=lambda0
    delta_lambda=10.0
    iter=0
    while ((abs(delta_lambda)>1e-9) && (iter<70))
        lambda0_qsn=CTa/(2*sqrt(mu^2+(lambda0_qs-muzs)^2))
        delta_lambda=lambda0_qsn-lambda0_qs
        lambda0_qs=lambda0_qs+0.5*delta_lambda
        iter=iter+1
    end

    # 计算迎角和质量流参数
    mutot=sqrt(mu^2+(muzs-lambda0_qs)^2)
    Vmpp=(mu^2+(lambda0_qs-muzs)*(2*lambda0_qs-muzs))/mutot

    # 风轴系到桨轴系的坐标转换
    cbetawind=muxs/mu
    sbetawind=muys/mu
    Ts2w=[1  0         0
          0  cbetawind sbetawind
          0 -sbetawind cbetawind]
    Tw2s=Ts2w'
    Tw2sdot=[0  0          0
             0 -sbetawind -cbetawind
             0 -cbetawind  sbetawind]*betawind_dot
    Ts2wdot=Tw2sdot'

    # 将入流状态、角速率、挥舞速率转换到风轴系
    lambda_w=Ts2w*[lambda0;lambda1s;lambda1c]
    omega_w=Ts2w*[rs;ps;qs]/OMEGA
    beta1c_dot=xr[NB+3]
    beta1s_dot=xr[NB+4]
    beta0_dot=xr[NB+1]
    betastar_w=Ts2w*[beta0_dot;beta1s_dot;beta1c_dot]/OMEGA

    # Pitt-Peters 模型，包含准稳态尾迹曲率效应 (基于 Zhao 2004)
    Xpp=xr[4*NB+6]
    kcpp=xr[4*NB+7]
    kspp=xr[4*NB+8]

    # 限制分母为正
    kcqs=(omega_w[3]-betastar_w[3])/max((lambda0_qs-muzs),0.001)
    ksqs=(omega_w[2]-betastar_w[2])/max((lambda0_qs-muzs),0.001)
    chi=atan2(mu,(lambda0_qs-muzs))
    Xqs=tan(0.5*chi)

    # 计算入流模型状态的导数
    Xpp_dot=OMEGA*(15.0*pi*Vmpp)/32.0*(Xqs-Xpp)
    kcpp_dot=OMEGA*(15.0*pi*mutot)/32.0*(kcqs-kcpp)
    kspp_dot=OMEGA*(15.0*pi*mutot)/32.0*(ksqs-kspp)

    # 定义 Pitt-Peters 模型的 Vpp 和 Lpp 矩阵
    Vpp=[mutot 0.0  0.0
         0.0   Vmpp 0.0
         0.0   0.0  Vmpp]
    Lpp=[0.5          0.0         -15*pi/64*Xpp
         0.0          2*(1+Xpp^2)  0
         15*pi/64*Xpp 0.0          2*(1-Xpp^2)]

    # 尾迹曲率修正 (当前假设为0)
    delL1=zeros(3,3)
    delL1[2,1]=0.5*KR*kspp
    delL1[3,1]=0.5*KR*kcpp
    delL1[1,2]=0.5*KR*kspp
    delL1[1,3]=0.5*KR*kcpp
    delL2=zeros(3,3)
    delL2[2,1]=0.75*KR*kspp*Xpp^2
    delL2[3,1]=-0.75*KR*kcpp*Xpp^2
    delL3=zeros(3,3)
    delL3[2,1]=1.25*mu*kcpp*Xpp*KR
    delL3[2,2]=(-2.5*kcpp*Xpp-1.5*mu*kspp*(1.0+1.5*Xpp^2))*KR
    delL3[2,3]=-2.5*kspp*Xpp*KR
    delL3[3,1]=1.25*mu*kspp*Xpp*KR
    delL3[3,2]=(-2.5*kspp*Xpp-1.5*mu*kcpp*(1.0-1.5*Xpp^2))*KR
    delL3[3,3]=-0.3*kcpp*Xpp*KR

    # 计算入流状态导数 (lambda_dot)
    ciao=OMEGA*[(75*pi)/128 0.0        0.0
                0.0         (45*pi)/16 0.0
                0.0         0.0        (45*pi)/16]
    ciao2=Ts2w*[CTa;-CLa;-CMa]-(inv(Lpp+delL1+delL2+delL3)*Vpp)*lambda_w
    lambda_dot_w=ciao*ciao2
    lambda_dot=Tw2s*lambda_dot_w+Tw2sdot*lambda_w
    # --- 7.8 状态导数和总力/力矩合成 ---

    # 构建主旋翼状态导数向量 xrdot
    xrdot=zeros(NRSTATES)
    xrdot[1:NB]=LIBC2MBC*beta_dot+LIBC2MBCdot*beta
    xrdot[NB+1:2*NB]=LIBC2MBC*beta_ddot+2*LIBC2MBCdot*beta_dot+LIBC2MBCddot*beta
    xrdot[2*NB+1:3*NB]=LIBC2MBC*zeta_dot+LIBC2MBCdot*zeta
    xrdot[3*NB+1:4*NB]=LIBC2MBC*zeta_ddot+2*LIBC2MBCdot*zeta_dot+LIBC2MBCddot*zeta
    xrdot[4*NB+1:4*NB+3]=lambda_dot
    xrdot[4*NB+4]=OMEGA
    xrdot[4*NB+5]=(1.0/TauDT)*(Fp_DT-xr[4*NB+5]) # 动态扭转力的滤波状态
    xrdot[4*NB+6]=Xpp_dot
    xrdot[4*NB+7]=kcpp_dot
    xrdot[4*NB+8]=kspp_dot

    # 将旋翼力和力矩转换到机体坐标系，并移动到飞机重心
    Fmrb=Tshaft'*[-Hforce;-Jforce;-Thrust];
    Mmrb=Tshaft'*[Lhub;Mhub;Qhub]+[0 -Zh Yh;Zh 0 -Xh;-Yh Xh 0]*Fmrb
    Fr=Fmrb
    Mr=Mmrb
    Qr=Qhub

    # ------------------------------------------------------------------------------
    #
    # ## 8. 运动方程求解
    #
    # ------------------------------------------------------------------------------

    # 合成作用在飞机上的总外力和外力矩
    Ftot=Fr+Ff+Fht+Fvt+Ftr
    Mtot=Mr+Mf+Mht+Mvt+Mtr
    Qreq=Qr+Qtr*GEARTR # 总需用功率

    # 调用运动方程模块，计算刚体运动的加速度
    xfdot=eqnmot(xf,Ftot,Mtot)

    # 组合成完整的飞机状态导数向量 xdot
    xdot=zeros(NSTATES)
    xdot[1:NFSTATES]=xfdot
    xdot[NFSTATES+1:NFSTATES+NRSTATES]=xrdot
    xdot[NSTATES]=xtrdot
    return xdot
end
