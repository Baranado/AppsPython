import math
def traslateSection(sigmac,fy,xgo,ygo,xso,yso,aso,kg,ks):
	xyp=centroidPlastic(sigmac,fy,xgo,ygo,xso,yso,aso,kg,ks)
	for i in range(kg):
		xg[i]=xgo[i]-xyp[0]
		yg[i]=ygo[i]-xyp[1]
	for i in range(ks):
		xs[i]=xso[i]-xyp[0]
		ys[i]=yso[i]-xyp[1]
	return [xg,yg,xs,ys,xyp[2]]
def centroidPlastic(sigmac,fy,xgo,ygo,xso,yso,aso,kg,ks):
	asxyg=momentPoly(xgo,ygo,kg)
	asxys=momentPunt(xso,yso,aso,ks)
	nyc=fy/sigmac+1
	at=asxyg[0]-nyc*asxys[0]
	xp=(asxyg[2]-nyc*asxys[2])/at
	yp=(asxyg[1]-nyc*asxys[1])/at
	nmin=0.8*phic*sigmac*at
	return [xp,yp,nmin]
def stressRebar(epsi,ncy,cond):
	if abs(epsi)<1:
		sigmas=epsi
	else:
		sigmas=math.copysign(1,epsi)
	if cond:
		sigmas+=ncy
	return sigmas
def factorSecurity(epsi,epsity):
	if epsi<=1:
		phik=phic
	elif epsi>=epsity:
		phik=phit
	else:
		phik=phic+(phit-phic)*(epsi-1)/(epsity-1)
	return phik
def momentPoly(x,y,k):
	a=[]
	Ao=0
	Sx=0
	Sy=0
	for i in range(k):
		a.append(x[i]*y[i+1]-x[i+1]*y[i])
		Ao+=a[i]
		Sx+=(y[i]+y[i+1])*a[i]
		Sy+=(x[i]+x[i+1])*a[i]
	return [Ao/2,Sx/6,Sy/6]
def areaPoly(x,y,k):
	Ao=0
	for i in range(k):
		Ao+=x[i]*y[i+1]-x[i+1]*y[i]
	return Ao/2
def momentPunt(x,y,a,k):
	Ao=0
	Sx=0
	Sy=0
	for i in range(k):
		Ao+=a[i]
		Sx+=y[i]*a[i]
		Sy+=x[i]*a[i]
	return [Ao,Sx,Sy]
if __name__ == '__main__':
	kg=4
	xgo=[0,9,26,7,0]
	ygo=[0,-3,19,24,0]
	ks=3
	xso=[7,18,10]
	yso=[3,16,18]
	aso=[3.14159265359,2.0106192983,2.0106192983]
	#print(momentPunt(xso,yso,aso,ks))
	#print(areaPoly(xgo,ygo,kg))
	print(momentPoly(xgo,ygo,kg))

