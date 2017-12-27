# TLVCsys, SysCond, DiaMat

TLVCsys:=proc(n::integer,u::symbol,v::symbol,r::symbol,s::symbol,d::symbol,e::symbol,
	a11::symbol:=a[11],a12::symbol:=a[12],a21::symbol:=a[21],a22::symbol:=a[22])
	local f,g,vars,M,Cond,i,j;
	description "Generate a two-species 2n-dimensional Lotka-Volterra competitive system with discrete diffusion";
	
	# System
	f:=seq(u[i]*(r[i]-a11*u[i]-a12*v[i])+add(d[i,j]*u[j],j=1..n),i=1..n);
	g:=seq(v[i]*(s[i]-a21*u[i]-a22*v[i])+add(e[i,j]*v[j],j=1..n),i=1..n);
	
	# Variables
	vars:=seq(u[i],i=1..n),seq(v[i],i=1..n);
	
	# Jacobian Matrix
	M:=VectorCalculus:-Jacobian([f,g],[vars]);
	
	# Equilibrium condition
	Cond:=seq(r[i]=-add(d[i,j]*u[j],j=1..n)/u[i]+a11*u[i]+a12*v[i],i=1..n),
		seq(s[i]=-add(e[i,j]*v[j],j=1..n)/v[i]+a21*u[i]+a22*v[i],i=1..n);
	
	# Substitue the equilibrium condition to the Jacobian Matrix and simplify
	simplify(subs(Cond,M));
end proc:

(*
 * # Example
 * TLVCsys(2, u, v, r, s, d, e);
 * TLVCsys(3, u, v, r, s, d, e);
 *)

SysCond:=proc(n::integer, d::symbol, e::symbol, cd::equation)
	description "Generate the symmetry condition and positive definite condition.";
	local i,j;
	seq(seq(d[j,i]=d[i,j],j=i+1..n),i=1..n-1),seq(seq(e[j,i]=e[i,j],j=i+1..n),i=1..n-1),cd;
end proc:

(*
 * # Example
 * SysCond(3, d, e, a[22] = epsilon+a[12]*a[21]/a[11]);
 *)

DiaMat:=proc(n::integer, u::symbol, v::symbol, a11::symbol:=a[11])
	description "Generate a diagonal matrix.";
	local i,j;
	LinearAlgebra:-DiagonalMatrix([seq(u[i],i=1..n),seq(a11*v[i],i=1..n)]);
end proc:

(*
 * # Example
 * DiaMat(5, u, v);
 *)

# Example
(*
n := 5;
M1 := TLVCsys(n, u, v, r, s, d, e):
M2 := simplify(subs(SysCond(n, d, e, a[22] = epsilon+a[12]*a[21]/a[11]), M1)):
M3 := DiaMat(n, u, v) . M2:

StringTools:-FormatTime("%Y-%m-%d,%X");
bt := time[real]():
F1 := [seq((-1)^i*LinearAlgebra:-Determinant(LinearAlgebra:-SubMatrix(M3, [1 .. i], [1 .. i])), i = 1 .. 2*n)]:
save F1, "BKSystem_F1.m";
et := time[real]()-bt;
StringTools:-FormatTime("%Y-%m-%d,%X");
*)

# ConstVars

ConstVars:=proc(n::integer, d::symbol, e::symbol, cd::symbol,
	a11::symbol:=a[11],a12::symbol:=a[12],a21::symbol:=a[21],a22::symbol:=a[22])
	description "Generate a list of variables.";
	local i,j;
	[cd, a11, a12, a21, a22, seq(seq(d[i,j],j=i+1..n),i=1..n-1), seq(seq(e[i,j],j=i+1..n),i=1..n-1)];
end proc:

(*
 * # Example
 * ConstVars(5, d, e, epsilon);
 *)

# FilterAllPositiveCoeff

HasNegativeItem:=proc(p::polynom)
	if min(coeffs(p))<0 then
		return(true);
	else
		return(false);
	fi;
end proc:

FilterAllPositiveCoeff:=proc(p::{polynom,list})
	description "Remove all positive coefficient polynomial.";
	local Temp;
	
	Temp:=`if`(type(p, polynom),[expand(p)],expand(p));

	[selectremove(HasNegativeItem,Temp)];
end proc:

(*
 * # Example
 * FilterAllPositiveCoeff([x-y, x+y, x^2, 1, -4, -y^2-x]);
 *)

# FilterEOPower

EOPower:=proc(p::polynom)
	local Temp, Tempodd, Tempeven, i;
	Temp:=sqrfree(p);
	if Temp[1]<0 then 
		return([[p],[]]);
	fi;
	Tempodd:=NULL;
	Tempeven:=NULL;
	for i in Temp[2] do
		if type(i[2],odd) then
			Tempodd:=Tempodd, i[1];
		else
			Tempeven:=Tempeven, i[1];
		fi;
	od;
	[[Tempodd],[Tempeven]];
end proc:

(*
 * # Example
 * EOPower(u[1]*u[3]^2*u[4]^2*u[5]^2*(u[1]*v[2]+u[2]*v[1])*(u[1]*v[2]-u[2]*v[1])^2);
 *)

FilterEOPower:=proc(p::{polynom,list})
	description "Select the odd power factors.";
	local Temp;

	Temp:=`if`(type(p, polynom), [p], p);
	
	Temp:=map(EOPower,Temp);
	
	if nops(Temp)=0 then
		return([[],[]]);
	fi;
	
	ListTools:-Transpose(Temp);
end proc:

(*
 * FilterEOPower(u[3]^2*u[4]^2*u[5]^2*(u[1]*v[2]+u[2]*v[1])*(u[1]*v[2]-u[2]*v[1])^2);
 * FilterEOPower(u[3]*u[4]*u[5]^2*(u[1]^3*v[2]^3-u[1]^2*u[2]*v[1]*v[2]^2
 * 	-u[1]*u[2]^2*v[1]^2*v[2]+u[2]^3*v[1]^3+u[3]^3*v[1]^3+u[3]^3*v[2]^3
 * 	+u[4]^3*v[1]^3+u[4]^3*v[2]^3));
 * FilterEOPower([u[3]^2*u[4]^2*u[5]^2*(u[1]*v[2]+u[2]*v[1])*(u[1]*v[2]-u[2]*v[1])^2, 
 * 	u[3]*u[4]*u[5]^2*(u[1]^3*v[2]^3-u[1]^2*u[2]*v[1]*v[2]^2-u[1]*u[2]^2*v[1]^2*v[2]
 * 	+u[2]^3*v[1]^3+u[3]^3*v[1]^3+u[3]^3*v[2]^3+u[4]^3*v[1]^3+u[4]^3*v[2]^3)]);
 *)

# ItemClassByNegaVarSet

HasSameVariables:=proc(a::polynom,b::polynom)
	if indets(a)=indets(b) then
		return(true);
	else
		return(false);
	fi;
end proc:

SelectRemoveMostVarItem:=proc(p::list)
	local PC,PF;
	uses ListTools;
	PC:=[Categorize((x,y)->numelems(indets(x))=numelems(indets(y)),p)];
	PF:=FindMaximalElement(PC,(x,y)->nops(x[1])<nops(y[1]),'position');
	PF[1],Flatten(subsop(PF[2]=NULL,PC));
end proc:

(*
 * # Example
 * SelectRemoveTheGreatestItem([u[1]^2*u[3]^2*u[4]*u[5]*v[1]*v[3]*v[4]^2*v[5]^2, 
u[1]^2*u[3]*u[4]^2*u[5]*v[1]*v[3]^2*v[4]*v[5]^2, u[1]*u[3]^2*u[4]*u[5]^2*v[1]^2*v[3]*v[4]^2*v[5], 
u[1]*u[3]*u[4]^2*u[5]^2*v[1]^2*v[3]^2*v[4]*v[5], u[1]^3*u[3]^3*v[4]^3*v[5]^3, 
u[1]^3*u[4]^3*v[3]^3*v[5]^3, u[3]^3*u[5]^3*v[1]^3*v[4]^3, u[4]^3*u[5]^3*v[1]^3*v[3]^3, 
-u[1]^3*u[3]^2*u[4]*v[3]*v[4]^2*v[5]^3, -u[1]^3*u[3]*u[4]^2*v[3]^2*v[4]*v[5]^3, 
-u[1]^2*u[3]^3*u[5]*v[1]*v[4]^3*v[5]^2, -u[1]^2*u[4]^3*u[5]*v[1]*v[3]^3*v[5]^2, 
-u[1]*u[3]^3*u[5]^2*v[1]^2*v[4]^3*v[5], -u[1]*u[4]^3*u[5]^2*v[1]^2*v[3]^3*v[5], 
-u[3]^2*u[4]*u[5]^3*v[1]^3*v[3]*v[4]^2, -u[3]*u[4]^2*u[5]^3*v[1]^3*v[3]^2*v[4]]);
 *)

MainAndAuxiVariables:=proc(p::list(polynom))
	local pFactor, MainVariables, pIndets, i, j;
	pFactor:=factor(convert(p,`+`));
	MainVariables:=seq(`if`(type(i, `+`),seq(indets(j), j in i),NULL),i in pFactor);
	pIndets:=indets(p) minus `union`(MainVariables);
	[MainVariables, `if`(pIndets={},NULL,pIndets)];
end proc:

(*
 * # Example
 * MainAndAuxiVariables([u[2]^2*u[3]^2*u[4]*u[5]*v[2]*v[3]*v[4]^2*v[5]^2, 
u[2]^2*u[3]*u[4]*u[5]^2*v[2]*v[3]^2*v[4]^2*v[5], 
u[2]*u[3]^2*u[4]^2*u[5]*v[2]^2*v[3]*v[4]*v[5]^2, 
u[2]*u[3]*u[4]^2*u[5]^2*v[2]^2*v[3]^2*v[4]*v[5]]);
 * MainAndAuxiVariables([u[2]^2*u[3]^2*u[4]*u[5]*v[2]*v[3]*v[4]^2*v[5]^2, 
u[2]^2*u[3]*u[4]*u[5]^2*v[2]*v[3]^2*v[4]^2*v[5]]);
 *)

IsInTheItemClass:=proc(p::monomial, MainAuxiVariableSets::list)
	local pIndets, i;
	pIndets:=indets(p);
	if not pIndets subset `union`(op(MainAuxiVariableSets)) then
		return(false);
	fi;
	for i in MainAuxiVariableSets do
		if i subset pIndets then
			pIndets:=pIndets minus i;
		fi;
		if pIndets={} then
			return(true);
		fi;
	od;
	false;
end proc:

(*
 * # Example
 * IsInTheItemClass(u[2]^3*u[3]^3*v[2]^3*v[5]^3, [{u[3], v[5]}, {u[5], v[3]}, {u[2], v[4]}, {u[4], v[2]}]);
 * IsInTheItemClass(u[2]^3*u[5]^3*v[2]^3*v[3]^3, [{u[3], v[5]}, {u[5], v[3]}, {u[2], v[4]}, {u[4], v[2]}]);
 * map(IsInTheItemClass, [op(u[2]^3*u[3]^3*v[2]^3*v[5]^3+u[2]^3*u[3]^3*v[3]^3*v[4]^3
+u[2]^3*u[5]^3*v[2]^3*v[3]^3+u[2]^3*u[5]^3*v[4]^3*v[5]^3+u[3]^3*u[4]^3*v[2]^3*v[3]^3
+u[3]^3*u[4]^3*v[4]^3*v[5]^3+u[4]^3*u[5]^3*v[2]^3*v[5]^3+u[4]^3*u[5]^3*v[3]^3*v[4]^3)], 
[{u[3], v[5]}, {u[5], v[3]}, {u[2], v[4]}, {u[4], v[2]}]);
 * map(IsInTheItemClass, [op(u[2]^3*u[3]^3*v[4]^3*v[5]^3-u[2]^3*u[3]^2*u[5]*v[3]*v[4]^3*v[5]^2
-u[2]^3*u[3]*u[5]^2*v[3]^2*v[4]^3*v[5]+u[2]^3*u[5]^3*v[3]^3*v[4]^3
-u[2]^2*u[3]^3*u[4]*v[2]*v[4]^2*v[5]^3+u[2]^2*u[3]^2*u[4]*u[5]*v[2]*v[3]*v[4]^2*v[5]^2
+u[2]^2*u[3]*u[4]*u[5]^2*v[2]*v[3]^2*v[4]^2*v[5]-u[2]^2*u[4]*u[5]^3*v[2]*v[3]^3*v[4]^2
-u[2]*u[3]^3*u[4]^2*v[2]^2*v[4]*v[5]^3+u[2]*u[3]^2*u[4]^2*u[5]*v[2]^2*v[3]*v[4]*v[5]^2
+u[2]*u[3]*u[4]^2*u[5]^2*v[2]^2*v[3]^2*v[4]*v[5]-u[2]*u[4]^2*u[5]^3*v[2]^2*v[3]^3*v[4]
+u[3]^3*u[4]^3*v[2]^3*v[5]^3-u[3]^2*u[4]^3*u[5]*v[2]^3*v[3]*v[5]^2
-u[3]*u[4]^3*u[5]^2*v[2]^3*v[3]^2*v[5]+u[4]^3*u[5]^3*v[2]^3*v[3]^3)], 
[{u[3], v[5]}, {u[5], v[3]}, {u[2], v[4]}, {u[4], v[2]}]);
 *)

ItemClassByGreat:=proc(p::polynom)
	description "Classify the item in a polynom by its greatest items' variables.";
	local Temp, GreatestItems, GreatestItemsClass, NormalItems, VariableSets,
		IsTheFirstNonNegativeItem, Out, i, j, iItem, iCoeff, jCoeff, tempout;

	if type(p,monomial) then
		return([[p]]);
	fi;

	Temp:=[op(expand(p))];

	GreatestItems, NormalItems:=SelectRemoveMostVarItem(Temp);
	GreatestItemsClass:=[ListTools:-Categorize(HasSameVariables,GreatestItems)];

	VariableSets:=map(MainAndAuxiVariables, GreatestItemsClass);
	
	Out:=GreatestItemsClass;
	
	IsTheFirstNonNegativeItem:=true;
	for i in NormalItems do
		iItem:=i;
		for j from 1 to nops(GreatestItemsClass) do
			if IsInTheItemClass(iItem, VariableSets[j]) then
				iCoeff:=abs(coeffs(iItem));
				jCoeff:=abs(coeffs(GreatestItemsClass[j][1]));
				if iCoeff>=jCoeff then
					tempout:=iItem*jCoeff/iCoeff;
				else
					tempout:=iItem;
				fi;
				iItem:=iItem-tempout;
				Out:=subsop(j=[op(Out[j]),tempout],Out);
			fi;
			if iItem=0 then
				break;
			fi;
		od;
		if IsTheFirstNonNegativeItem then
			Out:=[op(Out),[iItem]];
			IsTheFirstNonNegativeItem:=false;
		else
			Out:=subsop(-1=[op(Out[-1]),iItem],Out);
		fi;
	od;
	Out;
end proc:

(*
 * # Example
 * ItemClassByGreat(u[1]^3*u[2]^3*v[3]^3-u[1]^3*u[2]^2*u[3]*v[2]*v[3]^2-
 * 	u[1]^3*u[2]*u[3]^2*v[2]^2*v[3]+u[1]^3*u[3]^3*v[2]^3+u[2]^3*u[5]^3*v[1]^3+
 * 	u[2]^3*u[5]^3*v[3]^3-u[2]^2*u[3]*u[5]^3*v[2]*v[3]^2-u[2]*u[3]^2*u[5]^3*v[2]^2*v[3]+
 * 	u[3]^3*u[5]^3*v[1]^3+u[3]^3*u[5]^3*v[2]^3);
 * 
 * Test:=u[2]^3*u[3]^3*v[1]^3*v[4]^3+u[2]^3*u[3]^3*v[2]^3*v[4]^3
+u[2]^3*u[3]^3*v[3]^3*v[5]^3+2*u[2]^3*u[3]^3*v[4]^3*v[5]^3
-u[2]^3*u[3]^2*u[4]*v[1]^3*v[3]*v[4]^2-u[2]^3*u[3]^2*u[4]*v[2]^3*v[3]*v[4]^2
-u[2]^3*u[3]^2*u[4]*v[3]*v[4]^2*v[5]^3-u[2]^3*u[3]*u[4]^2*v[1]^3*v[3]^2*v[4]
-u[2]^3*u[3]*u[4]^2*v[2]^3*v[3]^2*v[4]-u[2]^3*u[3]*u[4]^2*v[3]^2*v[4]*v[5]^3
+u[2]^3*u[4]^3*v[1]^3*v[3]^3+u[2]^3*u[4]^3*v[2]^3*v[3]^3
+2*u[2]^3*u[4]^3*v[3]^3*v[5]^3+u[2]^3*u[4]^3*v[4]^3*v[5]^3
-u[2]^2*u[3]^3*u[5]*v[2]*v[3]^3*v[5]^2-u[2]^2*u[3]^3*u[5]*v[2]*v[4]^3*v[5]^2
-u[2]^2*u[4]^3*u[5]*v[2]*v[3]^3*v[5]^2-u[2]^2*u[4]^3*u[5]*v[2]*v[4]^3*v[5]^2
-u[2]*u[3]^3*u[5]^2*v[2]^2*v[3]^3*v[5]-u[2]*u[3]^3*u[5]^2*v[2]^2*v[4]^3*v[5]
-u[2]*u[4]^3*u[5]^2*v[2]^2*v[3]^3*v[5]-u[2]*u[4]^3*u[5]^2*v[2]^2*v[4]^3*v[5]
+u[3]^3*u[5]^3*v[1]^3*v[4]^3+u[3]^3*u[5]^3*v[2]^3*v[3]^3
+2*u[3]^3*u[5]^3*v[2]^3*v[4]^3+u[3]^3*u[5]^3*v[4]^3*v[5]^3
-u[3]^2*u[4]*u[5]^3*v[1]^3*v[3]*v[4]^2-u[3]^2*u[4]*u[5]^3*v[2]^3*v[3]*v[4]^2
-u[3]^2*u[4]*u[5]^3*v[3]*v[4]^2*v[5]^3-u[3]*u[4]^2*u[5]^3*v[1]^3*v[3]^2*v[4]
-u[3]*u[4]^2*u[5]^3*v[2]^3*v[3]^2*v[4]-u[3]*u[4]^2*u[5]^3*v[3]^2*v[4]*v[5]^3
+u[4]^3*u[5]^3*v[1]^3*v[3]^3+2*u[4]^3*u[5]^3*v[2]^3*v[3]^3
+u[4]^3*u[5]^3*v[2]^3*v[4]^3+u[4]^3*u[5]^3*v[3]^3*v[5]^3:
 *
 * Test:=(u[3]*v[4]+u[4]*v[3])*(u[3]*v[4]-u[4]*v[3])^2*((u[1]*v[5]+u[5]*v[1])
*(u[1]*v[5]-u[5]*v[1])^2+2*u[2]^3*v[1]^3+2*u[2]^3*v[5]^3);
 * 
 * Test := expand((x-y)^3+(t+z)^3);
 *)

# IsEquivalent

RememberPermute:=proc(n::integer)
	option remember;
	combinat:-permute(n);
end proc:

IsEquivalent:=proc(p::polynom, q::polynom)
	local vp,vq,vs,vr,co,i,dvp,dvq,ip,iq,idp,idq;

	if nops(p)<>nops(q) then
		return(false);
	fi;

	if degree(p)<>degree(q) then
		return(false);
	fi;
	
	vp:=convert(indets(p),list);
	vq:=convert(indets(q),list);
	
	if nops(vp)<>nops(vq) then
		return(false);
	fi;

	dvp:=sort(map[2](degree,p,vp));
	dvq:=sort(map[2](degree,q,vq));
	if dvp<>dvq then
		return(false);
	fi;

	ip:=[op(p)];
	iq:=[op(q)];
	idp:=sort(map(degree,ip));
	idq:=sort(map(degree,iq));
	if idp<>idq then
		return(false);
	fi;

	vs:=RememberPermute(nops(vq));

	for i in vs do
		vr:=map(op,i,vq);
		co:=zip(`=`,vr,vp);
		if Testzero(p-subs(co,q)) then
			return(true);
		fi;
	od;
	false;
end proc:

(*
 * # Example
 * IsEquivalent(a-b, c-d);
 * 
 *)

# PDSimplify

PDFilter:=proc(p::{polynom,list})
	local PDFilterTemp, PDF_FAPC_p, PDF_FEOP_p, PDF_MFAPC_ppp, PDF_MOP_pp, PDF_MCO_p;
	PDFilterTemp:=`if`(type(p, polynom), [p], p);

	PDF_FAPC_p:=FilterAllPositiveCoeff(PDFilterTemp)[1];

	PDF_FEOP_p:=FilterEOPower(PDF_FAPC_p)[1];

	PDF_MFAPC_ppp:=map(FilterAllPositiveCoeff, PDF_FEOP_p);
	
	PDF_MOP_pp:=map[2](op, 1, PDF_MFAPC_ppp);

	PDF_MCO_p:=map(convert, PDF_MOP_pp, `*`);

	FilterAllPositiveCoeff(PDF_MCO_p)[1];
end proc:

(*
 * # Example
 * PDFilter([u[3]^2*u[4]^2*u[5]^2*(u[1]*v[2]+u[2]*v[1])*(u[1]*v[2]-u[2]*v[1])^2, 
 u[3]*u[4]*u[5]^2*(u[1]^3*v[2]^3-u[1]^2*u[2]*v[1]*v[2]^2-u[1]*u[2]^2*v[1]^2*v[2]
 +u[2]^3*v[1]^3+u[3]^3*v[1]^3+u[3]^3*v[2]^3+u[4]^3*v[1]^3+u[4]^3*v[2]^3)]);
 *)

PDSimplify:=proc(p::{polynom,list})
	local PDSSimplifyTemp, PDS_PDF_p, PDS_MICN_ppp, PDS_MMCO_pp,
	PDS_MPDF_pp,PDS_MCO_p;
	
	PDSSimplifyTemp:=`if`(type(p, polynom), [p], p);

	PDS_PDF_p:=PDFilter(PDSSimplifyTemp);
	
	PDS_MICN_ppp:=map(ItemClassByGreat, PDS_PDF_p);

	PDS_MMCO_pp:=map[2](map,convert,PDS_MICN_ppp,`+`);

	PDS_MPDF_pp:=map(PDFilter,PDS_MMCO_pp);

	PDS_MCO_p:=map(convert, PDS_MPDF_pp, `+`);

	FilterAllPositiveCoeff(PDS_MCO_p)[1];
end proc:

# PD

PD:=proc(p::polynom, vars::list(name))
	local Exp, Cop, Sip1, Sip2, Ecp;

	# Expand polynomial p.
	Exp:=expand(p);
	printf("There are %d items contained in the input polynomial.\n", nops(Exp));

	# Get the coefficients of polynomial p about variables vars.
	Cop:=[coeffs(Exp,vars)];
	printf("There are %d coefficients contained in the input polynomial.\n", nops(Cop));

	# Simplify these coefficient polynomials.
	Sip1:=PDSimplify(Cop);
	printf("There are %d polynomials after first simplify.\n", nops(Sip1));

	# Simplify these coefficient polynomials.
	Sip2:=PDSimplify(Sip1);
	printf("There are %d polynomials after second simplify.\n", nops(Sip2));

	# Equivalence class;
	Ecp:=[ListTools:-Categorize(IsEquivalent, Sip2)];

	map[2](op,1,Ecp);
end proc:


IsDDLVSysGloballyStable:=proc(n::integer)
	local OutPrintA, OutPrintB, M1, M2, M3, F1, F2, i, BT, Varbs, OutCondition;
	uses LinearAlgebra;

	print(StringTools:-FormatTime("%Y-%m-%d,%X"));
	BT := time[real]():
	# Construct the system for judgment and the leading principal minors.
	
	M1 := TLVCsys(n, u, v, r, s, d, e):
	M2 := simplify(subs(SysCond(n, d, e, a[22] = epsilon+a[12]*a[21]/a[11]), M1)):
	M3 := DiaMat(n, u, v) . M2:
	F1 := [seq((-1)^i*Determinant(SubMatrix(M3, [1 .. i], [1 .. i])), i = 1 .. 2*n)]:
	
	printf("It took %f seconds to calculate the leading principal minors.\n", time[real]()-BT);
	BT:=time[real]();
	# Judge whether the leading principal minors alternate in sign.
	
	Varbs := ConstVars(n, d, e, epsilon);
	F2 := map(PD, F1, Varbs):

	printf("It took %f seconds to judge whether the leading principal minors alternate in sign.\n", time[real]()-BT);
	
	OutCondition:=map(nops,F2);

	OutPrintA:=cat(
	"The %d-patch two-species Lotka-Volterra competitive discrete diffusion system ",
	"is globally stable.\n"
	);

	OutPrintB:=cat(
	"The %d-patch two-species Lotka-Volterra competitive discrete diffusion system ",
	"is not globally stable.\n",
	"Because the %ath leading principal minors of the system don't alternate in sign.\n"
	);

	if convert(OutCondition,`+`)=0 then
		printf(OutPrintA, n);
		print(StringTools:-FormatTime("%Y-%m-%d,%X"));
		return(true);
	else
		printf(OutPrintB, n, [seq(`if`(OutConditon[i]=0,NULL,i),i=1..nops(OutCondition))]);
		print(StringTools:-FormatTime("%Y-%m-%d,%X"));
		return(false);
	fi;
end proc:

IsDDLVSysGloballyStable(4);


