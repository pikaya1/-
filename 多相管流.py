
import math
class beggs_brill:
    def __init__(self,derp,N,D,H,fw,Qo,p0,Rp,rouo,bubp):
        self.H=H             # m
        self.derp=derp/100   # Mpa/100m to Mpa/m
        self.N=N             # .
        self.D=D*0.0254      # in to m
        self.eh=H/N          # m
        self.tep=2.9/100     # du/100m to du/m
        self.Qo=Qo/86400     # m3/d to m3/s
        self.p0=p0           # Mpa
        self.fw=fw           # .
        self.Ql=self.Qo/(1-fw)
        self.Qw=self.Ql*fw
        self.Rp=Rp           # .
        self.T0=20+273.15    # K
        self.xrouo=rouo      # .gamao
        self.bubp=bubp       # Mpa
        self.h=0
        self.xroug=0.7
        self.roua=1.293      # kg/m3
        self.rouo=rouo*1000  # kg/m3
        self.pba=0
        self.tba=0
        self.nwp=p0
        self.derps=[]
        self.rouw=1000
        self.ax=0.845
        self.bx=0.5351
        self.cx=0.0173
        self.dx=2.96
        self.ex=0.305
        self.fx=-0.4473
        self.gx=0.0978
    def ctpt(self):
        self.h=self.h+self.eh
        pba=(2*self.nwp+self.eh*self.derp)/2
        tba=(2*self.T0+self.tep*self.h)/2
        self.tba=tba
        self.pba=pba

        return pba,tba
    def ctrs(self,t,p):
        #t/K  p/Mpa
        API=141.5/self.xrouo-131.5
        A=0.0125*API-0.00091*(1.8*(t-273.15)+32) #du and mpa
        Rs=self.xroug/5.615*pow((7.9688*p+1.4)*pow(10,A),1.2048)
        if p>self.bubp:
            Rs=self.xroug/5.615*pow((7.9688*self.bubp+1.4)*pow(10,A),1.2048)
        return Rs
    def ctbo(self,Rs,t):
        #Rs/. t/K
        F=5.6146*Rs*pow(self.xroug/self.xrouo,0.5)+1.25*(1.8*(t-273.15)+32)
        Bo=0.9759+0.00012*pow(F,1.2)
        return Bo
    def ctz(self,P,t):
        #t/K,P/Mpa
        #Cranmer方法+试算法计算天然气压缩因子
        Tc=-82.5+273.15
        Pc=4.6
        Tr=(t)/Tc #K
        Pr=P/Pc
        Zs=1.41
        for i in range(0, 100):
            pr=0.27*Pr/(Tr*Zs)
            Z=1.0+(0.3151-1.0467/Tr-0.5783/pow(Tr,3))*pr+(0.5353-0.6123/Tr)*pow(pr,2)+(0.6815*pow(pr,2)/pow(Tr,3))
            if abs(Z-Zs)<0.005:
                return Z
            else:
                Zs=Z
        return Z
    def ctniang(self,z):
            self.pba=4
            Mg=28.97*self.xroug
            T=1.8*(self.tba)
            a=(9.4+0.02*Mg)/(209+19*Mg+T)*pow(T,1.5)
            b=3.5+986/T+0.01*Mg
            c=2.4-0.2*b
            print("dd")
            print(c)
            print(b)
            roug=(self.xroug*self.roua*self.T0*self.pba)/(z*self.tba*self.p0)
            print(roug)
            roug=roug*0.01
            ug=a*math.exp(b*pow(roug,c))*pow(10,-4)
            return ug
    def ctniano(self):
        #None
        return 3.5    # mpas
    def ctbiaozh(self):
        #None
        return 0.018  # N/m
    def work(self):

        self.ctpt()
        #self.pba=4.3
        #self.tba=313
        Rs=self.ctrs(self.tba,self.pba)
        Bo=self.ctbo(Rs,self.tba)

        r1=(self.xrouo*1000+self.roua*Rs*self.xroug)/Bo
        z=self.ctz(self.pba,self.tba)

        roug1=self.xroug*self.roua*self.pba*self.T0/(z*self.tba*0.101)


        qg=0.101*self.tba*z*(self.Rp-Rs)*self.Qo/(self.pba*self.T0)
        ql=self.Qo*Bo
        Ap = math.pi*pow(self.D/2,2)
        #print(ql)
        vs1=ql/Ap
        vsg=qg/Ap
        vsw=self.Qw/Ap
        vm=vs1+vsg+vsw
        Gl=ql*r1
        Gw=self.rouw*self.Qw
        Gg=roug1*qg
        Gm=Gl+Gg+Gw
        El=(ql+self.Qw)/(ql+self.Qw+qg)
        Nfr=pow(vm,2)/(9.81*self.D)
        print(Nfr)
        print(El)
        niano=self.ctniano()
        niang=self.ctniang(z)
        nianw=1-(self.tba-273.15-20)*0.01

        nianm=El*(niano*r1+nianw*self.rouw)/(r1+self.rouw)+niang*(1-El)
        print(El)
        print(Nfr)
        Nvl=vs1*pow((r1*(1-self.fw)+self.rouw*self.fw)/(9.81*self.ctbiaozh()),0.25)
        print(Nvl)
        L1=316*pow(El,0.302)
        L2=92.52*pow(10,-5)*pow(El,-2.4684)
        L3=0.1*pow(El,-1.4516)
        L4=0.5*pow(El,-6.738)

        print(L1)
        print(L2)
        print(L3)
        print(L4)
        print(El)
        print(Nfr)
        if (El<0.01 and Nfr<L1) or (El>0.01 and Nfr<L2):
            print("分离流")
            self.ax=0.98
            self.bx=0.4846
            self.cx=0.0868
            self.dx = 0.011
            self.ex = -3.768
            self.fx = 3.539
            self.gx = -1.614
        if (El>=0.01 and Nfr<=L3 and Nfr>L2):
            print("过渡流")
            self.ax = 0.845
            self.bx = 0.5351
            self.cx = 0.0173
            self.dx = 2.96
            self.ex = 0.305
            self.fx = -0.4473
            self.gx = 0.0978
        if (El>=0.01 and El<0.4 and Nfr<L1 and Nfr>L3) or (El>=0.4 and Nfr<L4 and Nfr>L3):
            print("间歇流")
            self.ax = 0.845
            self.bx = 0.5351
            self.cx = 0.0173
            self.dx = 2.96
            self.ex = 0.305
            self.fx = -0.4473
            self.gx = 0.0978
        if (El<0.4 and Nfr>=L1 ) or (El>=0.4 and Nfr>L4 ):
            print("分散流")
            self.ax = 1.065
            self.bx = 0.5929
            self.cx = 0.0609
        Hl0=self.ax*pow(El,self.bx)/pow(Nfr,self.cx)
        C=(1-El)*math.log(self.dx*pow(El,self.ex)*pow(Nvl,self.fx)*pow(Nfr,self.gx))
        print(C)
        print(Hl0)
        pusai=1+0.3*C
        Hlst=Hl0*pusai
        y=El/pow(Hlst,2)

        S=math.log(y)/(-0.0523+3.182*math.log(y)-0.8725*pow(math.log(y),2)+0.01853*pow(math.log(y),4))
        Nre=self.D*vm*(r1*El+roug1*(1-El))/(((niano*r1+nianw*self.rouw)/(r1+self.rouw))*El+niang*(1-El))
        print(Nre)
        lanb=0.0056+0.5/pow(Nre,0.32)
        lan=lanb*math.exp(S)
        r2=r1*(1-self.fw)+self.rouw*self.fw
        dpdz=pow(10,-6)*( 9.81*(r2*Hlst+roug1*(1-Hlst))+lan*Gm*vm/(2*self.D*Ap)  )/(  1-((r2*Hlst+roug1*(1-Hlst))*vm*vsg)/(self.pba*pow(10,6))  )

        if abs(self.eh*(dpdz-self.derp))>0.02:
          self.h=self.h-self.eh

          print("000")
          print(self.eh*(dpdz-self.derp))
          self.derp = dpdz
          return 0
        else:
          self.nwp=self.nwp+self.derp*self.eh
          self.derps.append(self.derp)
          print("111")
          print(self.derp)
          print(self.eh * (dpdz - self.derp))
          return 1


    def run(self):
      print("fagfasasfaf")
      print(self.nwp)
      for j in range(0,self.N):
       for i in range(0,10):
         print("+1")
         x=self.work()
         if x==1:
             print("derps")
             print(self.derps)
             print(len(self.derps))
             print(self.h)
             break
      allp=self.p0
      for i in self.derps:
          allp=allp+i*self.eh
      print("配产油量为50，井口油压为0.5时的井底流压：")
      print(allp)





example=beggs_brill(1.15,20,3,1720,0.5,50,0.5,50,0.86,10.8)
example.run()
