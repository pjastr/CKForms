#############################################################################
##
#W  ckforms.gi                  CKForms Package                  
##
##  Installation file for functions of the CKForms package.
##


###############################################################################
InstallGlobalFunction( NonCompactDimension, function(arg)
local G, P, dimP;
    G:=arg[1];
    P:=CartanDecomposition( G ).P;
    dimP:= Dimension(P);
    return dimP;
end);

###############################################################################
InstallGlobalFunction( RealRank, function(arg)
local G, rr;
    G:=arg[1];
    rr:=Dimension(CartanSubspace(G));
    return rr;
end);

###############################################################################
InstallGlobalFunction( AHypRankForSimple, function(arg)
local G, id,out,name, namesplit;
    G:=arg[1];
    if IsomorphismOfRealSemisimpleLieAlgebras(G,RealFormById("A",1,0))<>false then
        out:=1;
    elif IsomorphismOfRealSemisimpleLieAlgebras(G,RealFormById("A",3,0))<>false then
        out:=2;
    elif IsomorphismOfRealSemisimpleLieAlgebras(G,RealFormById("B",2,0))<>false then
        out:=2;
    elif IsomorphismOfRealSemisimpleLieAlgebras(G,RealFormById("C",3,0))<>false then
        out:=3;
    elif IsomorphismOfRealSemisimpleLieAlgebras(G,RealFormById("C",4,0))<>false then
        out:=4;
    elif IsomorphismOfRealSemisimpleLieAlgebras(G,RealFormById("D",4,0))<>false then
        out:=4;
    else
        id := IdRealForm(G);
        if id[1]="A" and id[2]=2 and id[3]=3 then
            out:=1;
        elif id[1]="A" and id[2]=2 and id[3]=0 then
            out:=1;
        elif id[1]="A" and id[2]=3 and id[3]=5 then
            out:=2;
        elif id[1]="A" and id[2]=3 and id[3]=0 then
            out:=2;
        elif id[1]="A" and id[2]=4 and id[3]=4 then
            out:=2;
        elif id[1]="A" and id[2]=4 and id[3]=0 then
            out:=2;
        elif id[1]="A" and id[2]=5 and id[3]=6 then
            out:=3;
        elif id[1]="A" and id[2]=5 and id[3]=0 then
            out:=3;
        elif id[1]="A" and id[2]=5 and id[3]=5 then
            out:=1;
        elif id[1]="A" and id[2]=6 and id[3]=5 then
            out:=3;
        elif id[1]="A" and id[2]=6 and id[3]=0 then
            out:=3;
        elif id[1]="A" and id[2]=7 and id[3]=7 then
            out:=4;
        elif id[1]="A" and id[2]=7 and id[3]=0 then
            out:=4;
        elif id[1]="A" and id[2]=7 and id[3]=6 then
            out:=2;
        elif id[1]="A" and id[2]=8 and id[3]=6 then
            out:=4;
        elif id[1]="A" and id[2]=8 and id[3]=0 then
            out:=4;
        elif id[1]="D" and id[2]=5 and id[3]=7 then
            out:=4;
        elif id[1]="D" and id[2]=5 and id[3]=0 then
            out:=4;
        elif id[1]="D" and id[2]=7 and id[3]=0 then
            out:=6;
        elif id[1]="D" and id[2]=7 and id[3]=0 then
            out:=6;
        elif id[1]="E" and id[2]=6 and id[3]=2 then
            out:=4;
        elif id[1]="E" and id[2]=6 and id[3]=5 then
            out:=1;
        elif id[1]="E" and id[2]=6 and id[3]=0 then
            out:=4;
        else
            out:=RealRank(G);
        fi;
    fi;
    return out;
end);

###############################################################################
InstallGlobalFunction( GetIdFromName, function(arg)
local sg, out,types, ranks,ids,j,k,l;
    out:=[];
    types:=["A","B","C","D","E","F","G"];
    for j in [1..Length(types)] do
        if (types[j]="A") then
            ranks:=[1..8];
            for k in ranks do
                ids:=[0..NumberRealForms(types[j],k)];
                for l in ids do
                    if (corelg.namesimple([types[j],k,l])=arg[1]) then
                        out:=[types[j],k,l];
                        break;
                    fi;
                od;
            od;
        elif (types[j]="B") then
            ranks:=[2..8];
            for k in ranks do
                ids:=[0..NumberRealForms(types[j],k)];
                for l in ids do
                    if (corelg.namesimple([types[j],k,l])=arg[1]) then
                        out:=[types[j],k,l];
                        break;
                    fi;
                od;
            od;
        elif (types[j]="C") then
            ranks:=[3..8];
            for k in ranks do
                ids:=[0..NumberRealForms(types[j],k)];
                for l in ids do
                    if (corelg.namesimple([types[j],k,l])=arg[1]) then
                        out:=[types[j],k,l];
                        break;
                    fi;
                od;
            od;
        elif (types[j]="D") then
            ranks:=[4..8];
            for k in ranks do
                ids:=[0..NumberRealForms(types[j],k)];
                for l in ids do
                    if (corelg.namesimple([types[j],k,l])=arg[1]) then
                        out:=[types[j],k,l];
                        break;
                    fi;
                od;
            od;
        elif (types[j]="E") then
            ranks:=[6..8];
            for k in ranks do
                ids:=[0..NumberRealForms(types[j],k)];
                for l in ids do
                    if (corelg.namesimple([types[j],k,l])=arg[1]) then
                        out:=[types[j],k,l];
                        break;
                    fi;
                od;
            od;
        elif (types[j]="F") then
            ranks:=[4];
            for k in ranks do
                ids:=[0..NumberRealForms(types[j],k)];
                for l in ids do
                    if (corelg.namesimple([types[j],k,l])=arg[1]) then
                        out:=[types[j],k,l];
                        break;
                    fi;
                od;
            od;
        elif (types[j]="G") then
            ranks:=[2];
            for k in ranks do
                ids:=[0..NumberRealForms(types[j],k)];
                for l in ids do
                    if (corelg.namesimple([types[j],k,l])=arg[1]) then
                        out:=[types[j],k,l];
                        break;
                    fi;
                od;
            od;
        fi;
    od;
    return out;
end);

###############################################################################
InstallGlobalFunction( AHypRank, function(arg)
local G,out, name, namesplit,j, parts;
    G:=arg[1];
    name:=NameRealForm(G);
    namesplit:=SplitString(name,"+");
    if (Length(namesplit)=1) then
        out:=AHypRankForSimple(G);
    else
        out:=0;
        parts:=[1..Length(namesplit)];
        if namesplit[Length(namesplit)]=" a torus of  1 non-compact dimensions" then
            parts:=[1..Length(namesplit)-1];
            namesplit[Length(namesplit)-1]:=ReplacedString(namesplit[Length(namesplit)-1]," ","");
        elif namesplit[Length(namesplit)]=" a torus of 1 compact dimensions" then
            parts:=[1..Length(namesplit)-1];
            namesplit[Length(namesplit)-1]:=ReplacedString(namesplit[Length(namesplit)-1]," ","");
        fi;
        for j in parts do
            out:=out+AHypRankForSimple(RealFormById(GetIdFromName(namesplit[j])));
        od;
    fi;
    return out;
end);

###############################################################################
InstallGlobalFunction( CompactDimension, function(arg)
local G, K, dimK;
    G:=arg[1];
    K:=CartanDecomposition( G ).K;
    dimK:= Dimension(K);
    return dimK;
end);


###############################################################################
InstallGlobalFunction( Separate, function(arg)
local s, n, i,j, t;
    s:=arg[1];
    t:=[];
    n:=Length(s);
    for j in [0..n] do
        i:=3*j+1;
        if i<= n then
            if i+4<=n then
                if s[i]=s[i+3] and s[i+1]=s[i+4] then
                    if i+5>n and i=1 then
                        t:= [s{ [i..i+4] },""];
                    elif i+5>n and i>1 then
                        t:=[s{ [i..i+4] },s{ [1..i-2] }];
                    elif i+5<=n and i=1 then
                        t:=[s{ [i..i+4] },s{ [i+6..n] }];
                    else  
                        t:= [s{ [i..i+4] },Concatenation(s{ [1..i-2] },s{ [i+5..n] })];
                    fi;
                    break;
                fi;
            fi;
        fi;
    od;
    return t;
end);



###############################################################################
InstallGlobalFunction( NicerSemisimpleSubalgebras, function(arg)
local G, s, sub, List1,List2, ListSub;
    G:=arg[1];
    s:= LieAlgebraAndSubalgebras(G);
    sub:= s.subalgs;
    List1:= List( sub, SemiSimpleType );
    List2:=Unique(List1);
    ListSub:=List2;
    return ListSub;
end);

###############################################################################
InstallGlobalFunction( NumberOfSubalgebraClasses, function(arg)
local G,H, s, sub, List1, i, len,num;
    G:=arg[1];
    H:=arg[2];
    s:= LieAlgebraAndSubalgebras(G);
    sub:= s.subalgs;
    List1:= List( sub, SemiSimpleType );
    num:=0;
    len:=Length(List1);
    for i in [1..len] do
        if List1[i]=H then
            num:=num+1;
        fi;
    od;
    return num;
end);


###############################################################################
InstallGlobalFunction( PotentialSubalgebras, function(arg)
local G,H1,H2,H3,H4,H5,H6,H7,H8,GG,type, rank, id,idH, str,len,t1,r1,t2,r2,t3,r3,t4,r4,t5,r5,t6,r6,t7,r7,t8,r8,i,j,k,k1,k2,k3,k4,k5,k6,ListSub,Sub, rankRG, rankRH, dimpG, dimpH, dimkG, dimkH, s1, st1, sr1,s2,st2,sr2,s3,st3,sr3,s4,st4,sr4,s5,st5,sr5,s6,st6,sr6,s7,st7,sr7;
    type:=arg[1];
    rank:=arg[2];
    id:=arg[3];
    G:=RealFormById(type, rank, id);
    GG:=Concatenation(type,String(rank));
    rankRG:=RealRank(G);
    dimpG:=NonCompactDimension(G);
    dimkG:=CompactDimension(G);
    ListSub:=NicerSemisimpleSubalgebras(GG);
    len:=Length(ListSub);
    Sub:=[];
    for i in [1..len] do
        str:=ListSub[i];
        t1:=[str[1]];
        r1:=IntChar(str[2])-48;
        if Length(str)>3 then
            t2:=[str[4]];
            r2:=IntChar(str[5])-48;
        fi;
        if Length(str)>6 then
            t3:=[str[7]];
            r3:=IntChar(str[8])-48;
        fi;
        if Length(str)>9 then
            t4:=[str[10]];
            r4:=IntChar(str[11])-48;
        fi;
        if Length(str)>12 then
            t5:=[str[13]];
            r5:=IntChar(str[14])-48;
        fi;
        if Length(str)>15 then
            t6:=[str[16]];
            r6:=IntChar(str[17])-48;
        fi;
        if Length(str)>18 then
            t7:=[str[19]];
            r7:=IntChar(str[20])-48;
        fi;
        if Length(str)>21 then
            t8:=[str[22]];
            r8:=IntChar(str[23])-48;
        fi;
        if Length(str) =2 then
            for k in [1..NumberRealForms(t1,r1)] do
                H1:=RealFormById(t1,r1,k);
                rankRH:=RealRank(H1);
                dimpH:=NonCompactDimension(H1);
                dimkH:=CompactDimension(H1);
                if rankRG> rankRH  then
                    if dimpG>= dimpH then
                        if dimkG> dimkH then
                            idH:=[rankRH,dimpH,t1,r1,k];
                            Add(Sub,idH);
                        fi;
                    fi;
                fi;
            od;
        elif Length(str) =5 then
            for k in [1..NumberRealForms(t1,r1)] do
                for j in [1..NumberRealForms(t2,r2)] do
                    H1:=RealFormById(t1,r1,k);
                    H2:=RealFormById(t2,r2,j); 
                    rankRH:=RealRank(H1)+RealRank(H2);
                    dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2);
                    dimkH:=CompactDimension(H1)+CompactDimension(H2);               
                    if rankRG> rankRH  then
                        if dimpG>= dimpH then
                            if dimkG> dimkH then
                                idH:=[rankRH,dimpH,t1,r1,k,t2,r2,j];
                                Add(Sub,idH);
                            fi;
                        fi;
                    fi;
                od;
            od;
            if Length(Separate(str))=2 then
                s1:=Separate(str)[1];
                st1:=[s1[1]];
                sr1:=IntChar(s1[2])-48;
                H1:=RealFormById(st1,sr1,0);
                rankRH:=RealRank(H1);
                dimpH:=NonCompactDimension(H1);
                dimkH:=CompactDimension(H1);
                if rankRG> rankRH  then
                    if dimpG>= dimpH then
                        if dimkG> dimkH then
                            idH:=[rankRH,dimpH,st1,sr1,0];
                            Add(Sub,idH);
                        fi;
                    fi;
                fi; 
            fi;
        elif Length(str) =8 then
            for k in [1..NumberRealForms(t1,r1)] do
                for j in [1..NumberRealForms(t2,r2)] do
                    for k1 in [1..NumberRealForms(t3,r3)] do
                        H1:=RealFormById(t1,r1,k);
                        H2:=RealFormById(t2,r2,j);
                        H3:=RealFormById(t3,r3,k1);
                        rankRH:=RealRank(H1)+RealRank(H2)+RealRank(H3);
                        dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2)+NonCompactDimension(H3);
                        dimkH:=CompactDimension(H1)+CompactDimension(H2)+CompactDimension(H3);          
                        if rankRG> rankRH  then
                            if dimpG>= dimpH then
                                if dimkG> dimkH then
                                    idH:=[rankRH,dimpH,t1,r1,k,t2,r2,j,t3,r3,k1];
                                    Add(Sub,idH);
                                fi;
                            fi;
                        fi;
                    od;
                od;
            od;
            if Length(Separate(str))=2 then
                s1:=Separate(str)[1];
                st1:=[s1[1]];
                sr1:=IntChar(s1[2])-48;
                s2:=Separate(str)[2];
                st2:=[s2[1]];
                sr2:=IntChar(s2[2])-48;
                for j in [1..NumberRealForms(st2,sr2)] do
                    H1:=RealFormById(st1,sr1,0);
                    H2:=RealFormById(st2,sr2,j); 
                    rankRH:=RealRank(H1)+RealRank(H2);
                    dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2);
                    dimkH:=CompactDimension(H1)+CompactDimension(H2);   
                    if rankRG> rankRH  then
                        if dimpG>= dimpH then
                            if dimkG> dimkH then
                                idH:=[rankRH,dimpH,st1,sr1,0,st2,sr2,j];
                                Add(Sub,idH);
                            fi;
                        fi;
                    fi; 
                od;
            fi;
        elif Length(str) =11 then
            for k in [1..NumberRealForms(t1,r1)] do
                for j in [1..NumberRealForms(t2,r2)] do
                    for k1 in [1..NumberRealForms(t3,r3)] do
                        for k2 in [1..NumberRealForms(t4,r4)] do
                            H1:=RealFormById(t1,r1,k);
                            H2:=RealFormById(t2,r2,j);
                            H3:=RealFormById(t3,r3,k1);
                            H4:=RealFormById(t4,r4,k2);
                            rankRH:=RealRank(H1)+RealRank(H2)+RealRank(H3)+RealRank(H4);
                            dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2)+NonCompactDimension(H3)+NonCompactDimension(H4);
                            dimkH:=CompactDimension(H1)+CompactDimension(H2)+CompactDimension(H3)+CompactDimension(H4);          
                            if rankRG> rankRH  then
                                if dimpG>= dimpH then
                                    if dimkG> dimkH then
                                        idH:=[rankRH,dimpH,t1,r1,k,t2,r2,j,t3,r3,k1,t4,r4,k2];
                                        Add(Sub,idH);
                                    fi;
                                fi;
                            fi;
                        od;
                    od;
                od;
            od;
            if Length(Separate(str))=2 then
                s1:=Separate(str)[1];
                st1:=[s1[1]];
                sr1:=IntChar(s1[2])-48;
                s2:=Separate(str)[2];
                st2:=[s2[1]];
                sr2:=IntChar(s2[2])-48;
                s3:=Separate(str)[2];
                st3:=[s3[4]];
                sr3:=IntChar(s3[5])-48;
                for j in [1..NumberRealForms(st2,sr2)] do
                    for k1 in [1..NumberRealForms(st3,sr3)] do
                        H1:=RealFormById(st1,sr1,0);
                        H2:=RealFormById(st2,sr2,j); 
                        H3:=RealFormById(st3,sr3,k1);
                        rankRH:=RealRank(H1)+RealRank(H2)+RealRank(H3);
                        dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2)+NonCompactDimension(H3);
                        dimkH:=CompactDimension(H1)+CompactDimension(H2)+CompactDimension(H3);   
                        if rankRG> rankRH  then
                            if dimpG>= dimpH then
                                if dimkG> dimkH then
                                    idH:=[rankRH,dimpH,st1,sr1,0,st2,sr2,j,st3,sr3,k1];
                                    Add(Sub,idH);
                                fi;
                            fi;
                        fi; 
                    od;
                od;
            fi;
            if Length(Separate(Separate(str)[2]))=2 then
                s1:=Separate(str)[1];
                st1:=[s1[1]];
                sr1:=IntChar(s1[2])-48;
                s2:=Separate(str)[2];
                st2:=[s2[1]];
                sr2:=IntChar(s2[2])-48;
                H1:=RealFormById(st1,sr1,0);
                H2:=RealFormById(st2,sr2,0); 
                rankRH:=RealRank(H1)+RealRank(H2);
                dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2);
                dimkH:=CompactDimension(H1)+CompactDimension(H2);   
                if rankRG> rankRH  then
                    if dimpG>= dimpH then
                        if dimkG> dimkH then
                            idH:=[rankRH,dimpH,st1,sr1,0,st2,sr2,0];
                            Add(Sub,idH);
                        fi;
                    fi;
                fi; 
            fi;
        elif Length(str) =14 then
            for k in [1..NumberRealForms(t1,r1)] do
                for j in [1..NumberRealForms(t2,r2)] do
                    for k1 in [1..NumberRealForms(t3,r3)] do
                        for k2 in [1..NumberRealForms(t4,r4)] do
                            for k3 in [1..NumberRealForms(t5,r5)] do
                                H1:=RealFormById(t1,r1,k);
                                H2:=RealFormById(t2,r2,j);
                                H3:=RealFormById(t3,r3,k1);
                                H4:=RealFormById(t4,r4,k2);
                                H5:=RealFormById(t5,r5,k3);
                                rankRH:=RealRank(H1)+RealRank(H2)+RealRank(H3)+RealRank(H4)+RealRank(H5);
                                dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2)+NonCompactDimension(H3)+NonCompactDimension(H4)+NonCompactDimension(H5);
                                dimkH:=CompactDimension(H1)+CompactDimension(H2)+CompactDimension(H3)+CompactDimension(H4)+CompactDimension(H5);          
                                if rankRG> rankRH  then
                                    if dimpG>= dimpH then
                                        if dimkG> dimkH then
                                            idH:=[rankRH,dimpH,t1,r1,k,t2,r2,j,t3,r3,k1,t4,r4,k2,t5,r5,k3];
                                            Add(Sub,idH);
                                        fi;
                                    fi;
                                fi;
                            od;
                        od;
                    od;
                od;
            od;
            if Length(Separate(str))=2 then
                s1:=Separate(str)[1];
                st1:=[s1[1]];
                sr1:=IntChar(s1[2])-48;
                s2:=Separate(str)[2];
                st2:=[s2[1]];
                sr2:=IntChar(s2[2])-48;
                s3:=Separate(str)[2];
                st3:=[s3[4]];
                sr3:=IntChar(s3[5])-48;
                s4:=Separate(str)[2];
                st4:=[s4[7]];
                sr4:=IntChar(s4[8])-48;
                for j in [1..NumberRealForms(st2,sr2)] do
                    for k1 in [1..NumberRealForms(st3,sr3)] do
                        for k2 in [1..NumberRealForms(st4,sr4)] do
                            H1:=RealFormById(st1,sr1,0);
                            H2:=RealFormById(st2,sr2,j); 
                            H3:=RealFormById(st3,sr3,k1);
                            H4:=RealFormById(st4,sr4,k2);
                            rankRH:=RealRank(H1)+RealRank(H2)+RealRank(H3)+RealRank(H4);
                            dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2)+NonCompactDimension(H3)+NonCompactDimension(H4);
                            dimkH:=CompactDimension(H1)+CompactDimension(H2)+CompactDimension(H3)+CompactDimension(H4);   
                            if rankRG> rankRH  then
                                if dimpG>= dimpH then
                                    if dimkG> dimkH then
                                        idH:=[rankRH,dimpH,st1,sr1,0,st2,sr2,j,st3,sr3,k1,st4,sr4,k2];
                                        Add(Sub,idH);
                                    fi;
                                fi;
                            fi; 
                        od;
                    od;
                od;
            fi;
            if Length(Separate(Separate(str)[2]))=2 then
                s1:=Separate(str)[1];
                st1:=[s1[1]];
                sr1:=IntChar(s1[2])-48;
                s2:=Separate(Separate(str)[2])[1];
                st2:=[s2[1]];
                sr2:=IntChar(s2[2])-48;
                s3:=Separate(Separate(str)[2])[2];
                st3:=[s3[1]];
                sr3:=IntChar(s3[2])-48;
                for k1 in [1..NumberRealForms(st3,sr3)] do
                    H1:=RealFormById(st1,sr1,0);
                    H2:=RealFormById(st2,sr2,0);
                    H3:=RealFormById(st3,sr3,k1);
                    rankRH:=RealRank(H1)+RealRank(H2)+RealRank(H3);
                    dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2)+NonCompactDimension(H3);
                    dimkH:=CompactDimension(H1)+CompactDimension(H2)+CompactDimension(H3);   
                    if rankRG> rankRH  then
                        if dimpG>= dimpH then
                            if dimkG> dimkH then
                                idH:=[rankRH,dimpH,st1,sr1,0,st2,sr2,0,st3,sr3,k1];
                                Add(Sub,idH);
                            fi;
                        fi;
                    fi; 
                od;
            fi;
        elif Length(str) =17 then
            for k in [1..NumberRealForms(t1,r1)] do
                for j in [1..NumberRealForms(t2,r2)] do
                    for k1 in [1..NumberRealForms(t3,r3)] do
                        for k2 in [1..NumberRealForms(t4,r4)] do
                            for k3 in [1..NumberRealForms(t5,r5)] do
                                for k4 in [1..NumberRealForms(t6,r6)] do
                                    H1:=RealFormById(t1,r1,k);
                                    H2:=RealFormById(t2,r2,j);
                                    H3:=RealFormById(t3,r3,k1);
                                    H4:=RealFormById(t4,r4,k2);
                                    H5:=RealFormById(t5,r5,k3);
                                    H6:=RealFormById(t6,r6,k4);
                                    rankRH:=RealRank(H1)+RealRank(H2)+RealRank(H3)+RealRank(H4)+RealRank(H5)+RealRank(H6);
                                    dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2)+NonCompactDimension(H3)+NonCompactDimension(H4)+NonCompactDimension(H5)+NonCompactDimension(H6);
                                    dimkH:=CompactDimension(H1)+CompactDimension(H2)+CompactDimension(H3)+CompactDimension(H4)+CompactDimension(H5)+CompactDimension(H6);          
                                    if rankRG> rankRH  then
                                        if dimpG>= dimpH then
                                            if dimkG> dimkH then
                                                idH:=[rankRH,dimpH,t1,r1,k,t2,r2,j,t3,r3,k1,t4,r4,k2,t5,r5,k3,t6,r6,k4];
                                                Add(Sub,idH);
                                            fi;
                                        fi;
                                    fi;
                                od;
                            od;
                        od;
                    od;
                od;
            od;
            if Length(Separate(str))=2 then
                s1:=Separate(str)[1];
                st1:=[s1[1]];
                sr1:=IntChar(s1[2])-48;
                s2:=Separate(str)[2];
                st2:=[s2[1]];
                sr2:=IntChar(s2[2])-48;
                s3:=Separate(str)[2];
                st3:=[s3[4]];
                sr3:=IntChar(s3[5])-48;
                s4:=Separate(str)[2];
                st4:=[s4[7]];
                sr4:=IntChar(s4[8])-48;
                s5:=Separate(str)[2];
                st5:=[s5[10]];
                sr5:=IntChar(s5[11])-48;
                for j in [1..NumberRealForms(st2,sr2)] do
                    for k1 in [1..NumberRealForms(st3,sr3)] do
                        for k2 in [1..NumberRealForms(st4,sr4)] do
                            for k3 in [1..NumberRealForms(st5,sr5)] do
                                H1:=RealFormById(st1,sr1,0);
                                H2:=RealFormById(st2,sr2,j); 
                                H3:=RealFormById(st3,sr3,k1);
                                H4:=RealFormById(st4,sr4,k2);
                                H5:=RealFormById(st5,sr5,k3);
                                rankRH:=RealRank(H1)+RealRank(H2)+RealRank(H3)+RealRank(H4)+RealRank(H5);
                                dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2)+NonCompactDimension(H3)+NonCompactDimension(H4)+NonCompactDimension(H5);
                                dimkH:=CompactDimension(H1)+CompactDimension(H2)+CompactDimension(H3)+CompactDimension(H4)+CompactDimension(H5);   
                                if rankRG> rankRH  then
                                    if dimpG>= dimpH then
                                        if dimkG> dimkH then
                                            idH:=[rankRH,dimpH,st1,sr1,0,st2,sr2,j,st3,sr3,k1,st4,sr4,k2,st5,sr5,k3];
                                            Add(Sub,idH);
                                        fi;
                                    fi;
                                fi; 
                            od;
                        od;
                    od;
                od;
            fi;
            if Length(Separate(Separate(str)[2]))=2 then
                s1:=Separate(str)[1];
                st1:=[s1[1]];
                sr1:=IntChar(s1[2])-48;
                s2:=Separate(Separate(str)[2])[1];
                st2:=[s2[1]];
                sr2:=IntChar(s2[2])-48;
                s3:=Separate(Separate(str)[2])[2];
                st3:=[s3[1]];
                sr3:=IntChar(s3[2])-48;
                s4:=Separate(Separate(str)[2])[2];
                st4:=[s4[4]];
                sr4:=IntChar(s4[5])-48;
                for k1 in [1..NumberRealForms(st3,sr3)] do
                    for k2 in [1..NumberRealForms(st4,sr4)] do
                        H1:=RealFormById(st1,sr1,0);
                        H2:=RealFormById(st2,sr2,0);
                        H3:=RealFormById(st3,sr3,k1);
                        H4:=RealFormById(st4,sr4,k2);
                        rankRH:=RealRank(H1)+RealRank(H2)+RealRank(H3)+RealRank(H4);
                        dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2)+NonCompactDimension(H3)+NonCompactDimension(H4);
                        dimkH:=CompactDimension(H1)+CompactDimension(H2)+CompactDimension(H3)+CompactDimension(H4);   
                        if rankRG> rankRH  then
                            if dimpG>= dimpH then
                                if dimkG> dimkH then
                                    idH:=[rankRH,dimpH,st1,sr1,0,st2,sr2,0,st3,sr3,k1,st4,sr4,k2];
                                    Add(Sub,idH);
                                fi;
                            fi;
                        fi; 
                    od;
                od;
            fi;
            if Length(Separate(Separate(Separate(str)[2])[2]))=2 then
                s1:=Separate(str)[1];
                st1:=[s1[1]];
                sr1:=IntChar(s1[2])-48;
                s2:=Separate(Separate(str)[2])[1];
                st2:=[s2[1]];
                sr2:=IntChar(s2[2])-48;
                s3:=Separate(Separate(str)[2])[2];
                st3:=[s3[1]];
                sr3:=IntChar(s3[2])-48;
                H1:=RealFormById(st1,sr1,0);
                H2:=RealFormById(st2,sr2,0);
                H3:=RealFormById(st3,sr3,0);
                rankRH:=RealRank(H1)+RealRank(H2)+RealRank(H3);
                dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2)+NonCompactDimension(H3);
                dimkH:=CompactDimension(H1)+CompactDimension(H2)+CompactDimension(H3);   
                if rankRG> rankRH  then
                    if dimpG>= dimpH then
                        if dimkG> dimkH then
                            idH:=[rankRH,dimpH,st1,sr1,0,st2,sr2,0,st3,sr3,0];
                            Add(Sub,idH);
                        fi;
                    fi;
                fi; 
            fi;
        elif Length(str) =20 then
            for k in [1..NumberRealForms(t1,r1)] do
                for j in [1..NumberRealForms(t2,r2)] do
                    for k1 in [1..NumberRealForms(t3,r3)] do
                        for k2 in [1..NumberRealForms(t4,r4)] do
                            for k3 in [1..NumberRealForms(t5,r5)] do
                                for k4 in [1..NumberRealForms(t6,r6)] do
                                    for k5 in [1..NumberRealForms(t7,r7)] do
                                        H1:=RealFormById(t1,r1,k);
                                        H2:=RealFormById(t2,r2,j);
                                        H3:=RealFormById(t3,r3,k1);
                                        H4:=RealFormById(t4,r4,k2);
                                        H5:=RealFormById(t5,r5,k3);
                                        H6:=RealFormById(t6,r6,k4);
                                        H7:=RealFormById(t7,r7,k5);
                                        rankRH:=RealRank(H1)+RealRank(H2)+RealRank(H3)+RealRank(H4)+RealRank(H5)+RealRank(H6)+RealRank(H7);
                                        dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2)+NonCompactDimension(H3)+NonCompactDimension(H4)+NonCompactDimension(H5)+NonCompactDimension(H6)+NonCompactDimension(H7);
                                        dimkH:=CompactDimension(H1)+CompactDimension(H2)+CompactDimension(H3)+CompactDimension(H4)+CompactDimension(H5)+CompactDimension(H6)+CompactDimension(H7);          
                                        if rankRG> rankRH  then
                                            if dimpG>= dimpH then
                                                if dimkG> dimkH then
                                                    idH:=[rankRH,dimpH,t1,r1,k,t2,r2,j,t3,r3,k1,t4,r4,k2,t5,r5,k3,t6,r6,k4,t7,r7,k5];
                                                    Add(Sub,idH);
                                                fi;
                                            fi;
                                        fi;
                                    od;
                                od;
                            od;
                        od;
                    od;
                od;
            od;
            if Length(Separate(str))=2 then
                s1:=Separate(str)[1];
                st1:=[s1[1]];
                sr1:=IntChar(s1[2])-48;
                s2:=Separate(str)[2];
                st2:=[s2[1]];
                sr2:=IntChar(s2[2])-48;
                s3:=Separate(str)[2];
                st3:=[s3[4]];
                sr3:=IntChar(s3[5])-48;
                s4:=Separate(str)[2];
                st4:=[s4[7]];
                sr4:=IntChar(s4[8])-48;
                s5:=Separate(str)[2];
                st5:=[s5[10]];
                sr5:=IntChar(s5[11])-48;
                s6:=Separate(str)[2];
                st6:=[s6[13]];
                sr6:=IntChar(s6[14])-48;
                for j in [1..NumberRealForms(st2,sr2)] do
                    for k1 in [1..NumberRealForms(st3,sr3)] do
                        for k2 in [1..NumberRealForms(st4,sr4)] do
                            for k3 in [1..NumberRealForms(st5,sr5)] do
                                for k4 in [1..NumberRealForms(st6,sr6)] do
                                    H1:=RealFormById(st1,sr1,0);
                                    H2:=RealFormById(st2,sr2,j); 
                                    H3:=RealFormById(st3,sr3,k1);
                                    H4:=RealFormById(st4,sr4,k2);
                                    H5:=RealFormById(st5,sr5,k3);
                                    H6:=RealFormById(st6,sr6,k4);
                                    rankRH:=RealRank(H1)+RealRank(H2)+RealRank(H3)+RealRank(H4)+RealRank(H5)+RealRank(H6);
                                    dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2)+NonCompactDimension(H3)+NonCompactDimension(H4)+NonCompactDimension(H5)+NonCompactDimension(H6);
                                    dimkH:=CompactDimension(H1)+CompactDimension(H2)+CompactDimension(H3)+CompactDimension(H4)+CompactDimension(H5)+CompactDimension(H6);   
                                    if rankRG> rankRH  then
                                        if dimpG>= dimpH then
                                            if dimkG> dimkH then
                                                idH:=[rankRH,dimpH,st1,sr1,0,st2,sr2,j,st3,sr3,k1,st4,sr4,k2,st5,sr5,k3,st6,sr6,k4];
                                                Add(Sub,idH);
                                            fi;
                                        fi;
                                    fi; 
                                od;
                            od;
                        od;
                    od;
                od;
            fi;
            if Length(Separate(Separate(str)[2]))=2 then
                s1:=Separate(str)[1];
                st1:=[s1[1]];
                sr1:=IntChar(s1[2])-48;
                s2:=Separate(Separate(str)[2])[1];
                st2:=[s2[1]];
                sr2:=IntChar(s2[2])-48;
                s3:=Separate(Separate(str)[2])[2];
                st3:=[s3[1]];
                sr3:=IntChar(s3[2])-48;
                s4:=Separate(Separate(str)[2])[2];
                st4:=[s4[4]];
                sr4:=IntChar(s4[5])-48;
                s5:=Separate(Separate(str)[2])[2];
                st5:=[s5[7]];
                sr5:=IntChar(s5[8])-48;
                for k1 in [1..NumberRealForms(st3,sr3)] do
                    for k2 in [1..NumberRealForms(st4,sr4)] do
                        for k3 in [1..NumberRealForms(st5,sr5)] do
                            H1:=RealFormById(st1,sr1,0);
                            H2:=RealFormById(st2,sr2,0);
                            H3:=RealFormById(st3,sr3,k1);
                            H4:=RealFormById(st4,sr4,k2);
                            H5:=RealFormById(st5,sr5,k3);
                            rankRH:=RealRank(H1)+RealRank(H2)+RealRank(H3)+RealRank(H4)+RealRank(H5);
                            dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2)+NonCompactDimension(H3)+NonCompactDimension(H4)+NonCompactDimension(H5);
                            dimkH:=CompactDimension(H1)+CompactDimension(H2)+CompactDimension(H3)+CompactDimension(H4)+CompactDimension(H5);   
                            if rankRG> rankRH  then
                                if dimpG>= dimpH then
                                    if dimkG> dimkH then
                                        idH:=[rankRH,dimpH,st1,sr1,0,st2,sr2,0,st3,sr3,k1,st4,sr4,k2,st5,sr5,k3];
                                        Add(Sub,idH);
                                    fi;
                                fi;
                            fi; 
                        od;
                    od;
                od;
            fi;
            if Length(Separate(Separate(Separate(str)[2])[2]))=2 then
                s1:=Separate(str)[1];
                st1:=[s1[1]];
                sr1:=IntChar(s1[2])-48;
                s2:=Separate(Separate(str)[2])[1];
                st2:=[s2[1]];
                sr2:=IntChar(s2[2])-48;
                s3:=Separate(Separate(Separate(str)[2])[2])[1];
                st3:=[s3[1]];
                sr3:=IntChar(s3[2])-48;
                s4:=Separate(Separate(Separate(str)[2])[2])[2];
                st4:=[s4[1]];
                sr4:=IntChar(s4[2])-48;
                for k2 in [1..NumberRealForms(st4,sr4)] do
                    H1:=RealFormById(st1,sr1,0);
                    H2:=RealFormById(st2,sr2,0);
                    H3:=RealFormById(st3,sr3,0);
                    H4:=RealFormById(st4,sr4,k2);
                    rankRH:=RealRank(H1)+RealRank(H2)+RealRank(H3)+RealRank(H4);
                    dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2)+NonCompactDimension(H3)+NonCompactDimension(H4);
                    dimkH:=CompactDimension(H1)+CompactDimension(H2)+CompactDimension(H3)+CompactDimension(H4);   
                    if rankRG> rankRH  then
                        if dimpG>= dimpH then
                            if dimkG> dimkH then
                                idH:=[rankRH,dimpH,st1,sr1,0,st2,sr2,0,st3,sr3,0,st4,sr4,k2];
                                Add(Sub,idH);
                            fi;
                        fi;
                    fi; 
                od;
            fi;
        elif Length(str) =23 then
            for k in [1..NumberRealForms(t1,r1)] do
                for j in [1..NumberRealForms(t2,r2)] do
                    for k1 in [1..NumberRealForms(t3,r3)] do
                        for k2 in [1..NumberRealForms(t4,r4)] do
                            for k3 in [1..NumberRealForms(t5,r5)] do
                                for k4 in [1..NumberRealForms(t6,r6)] do
                                    for k5 in [1..NumberRealForms(t7,r7)] do
                                        for k6 in [1..NumberRealForms(t8,r8)] do
                                            H1:=RealFormById(t1,r1,k);
                                            H2:=RealFormById(t2,r2,j);
                                            H3:=RealFormById(t3,r3,k1);
                                            H4:=RealFormById(t4,r4,k2);
                                            H5:=RealFormById(t5,r5,k3);
                                            H6:=RealFormById(t6,r6,k4);
                                            H7:=RealFormById(t7,r7,k5);
                                            H8:=RealFormById(t8,r8,k6);
                                            rankRH:=RealRank(H1)+RealRank(H2)+RealRank(H3)+RealRank(H4)+RealRank(H5)+RealRank(H6)+RealRank(H7)+RealRank(H8);
                                            dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2)+NonCompactDimension(H3)+NonCompactDimension(H4)+NonCompactDimension(H5)+NonCompactDimension(H6)+NonCompactDimension(H7)+NonCompactDimension(H8);
                                            dimkH:=CompactDimension(H1)+CompactDimension(H2)+CompactDimension(H3)+CompactDimension(H4)+CompactDimension(H5)+CompactDimension(H6)+CompactDimension(H7)+CompactDimension(H8);   
                                            if rankRG> rankRH  then
                                                if dimpG>= dimpH then
                                                    if dimkG> dimkH then
                                                        idH:=[rankRH,dimpH,t1,r1,k,t2,r2,j,t3,r3,k1,t4,r4,k2,t5,r5,k3,t6,r6,k4,t7,r7,k5,t8,r8,k6];
                                                        Add(Sub,idH);
                                                    fi;
                                                fi;
                                            fi;
                                        od;
                                    od;
                                od;
                            od;
                        od;
                    od;
                od;
            od;
            if Length(Separate(str))=2 then
                s1:=Separate(str)[1];
                st1:=[s1[1]];
                sr1:=IntChar(s1[2])-48;
                s2:=Separate(str)[2];
                st2:=[s2[1]];
                sr2:=IntChar(s2[2])-48;
                s3:=Separate(str)[2];
                st3:=[s3[4]];
                sr3:=IntChar(s3[5])-48;
                s4:=Separate(str)[2];
                st4:=[s4[7]];
                sr4:=IntChar(s4[8])-48;
                s5:=Separate(str)[2];
                st5:=[s5[10]];
                sr5:=IntChar(s5[11])-48;
                s6:=Separate(str)[2];
                st6:=[s6[13]];
                sr6:=IntChar(s6[14])-48;
                s7:=Separate(str)[2];
                st7:=[s7[16]];
                sr7:=IntChar(s7[17])-48;
                for j in [1..NumberRealForms(st2,sr2)] do
                    for k1 in [1..NumberRealForms(st3,sr3)] do
                        for k2 in [1..NumberRealForms(st4,sr4)] do
                            for k3 in [1..NumberRealForms(st5,sr5)] do
                                for k4 in [1..NumberRealForms(st6,sr6)] do
                                    for k5 in [1..NumberRealForms(st7,sr7)] do
                                        H1:=RealFormById(st1,sr1,0);
                                        H2:=RealFormById(st2,sr2,j); 
                                        H3:=RealFormById(st3,sr3,k1);
                                        H4:=RealFormById(st4,sr4,k2);
                                        H5:=RealFormById(st5,sr5,k3);
                                        H6:=RealFormById(st6,sr6,k4);
                                        H7:=RealFormById(st7,sr7,k5);
                                        rankRH:=RealRank(H1)+RealRank(H2)+RealRank(H3)+RealRank(H4)+RealRank(H5)+RealRank(H6)+RealRank(H7);
                                        dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2)+NonCompactDimension(H3)+NonCompactDimension(H4)+NonCompactDimension(H5)+NonCompactDimension(H6)+NonCompactDimension(H7);
                                        dimkH:=CompactDimension(H1)+CompactDimension(H2)+CompactDimension(H3)+CompactDimension(H4)+CompactDimension(H5)+CompactDimension(H6)+CompactDimension(H7);   
                                        if rankRG> rankRH  then
                                            if dimpG>= dimpH then
                                                if dimkG> dimkH then
                                                    idH:=[rankRH,dimpH,st1,sr1,0,st2,sr2,j,st3,sr3,k1,st4,sr4,k2,st5,sr5,k3,st6,sr6,k4,st7,sr7,k5];
                                                    Add(Sub,idH);
                                                fi;
                                            fi;
                                        fi; 
                                    od;
                                od;
                            od;
                        od;
                    od;
                od;
            fi;
            if Length(Separate(Separate(str)[2]))=2 then
                s1:=Separate(str)[1];
                st1:=[s1[1]];
                sr1:=IntChar(s1[2])-48;
                s2:=Separate(Separate(str)[2])[1];
                st2:=[s2[1]];
                sr2:=IntChar(s2[2])-48;
                s3:=Separate(Separate(str)[2])[2];
                st3:=[s3[1]];
                sr3:=IntChar(s3[2])-48;
                s4:=Separate(Separate(str)[2])[2];
                st4:=[s4[4]];
                sr4:=IntChar(s4[5])-48;
                s5:=Separate(Separate(str)[2])[2];
                st5:=[s5[7]];
                sr5:=IntChar(s5[8])-48;
                s6:=Separate(Separate(str)[2])[2];
                st6:=[s6[10]];
                sr6:=IntChar(s6[11])-48;
                for k1 in [1..NumberRealForms(st3,sr3)] do
                    for k2 in [1..NumberRealForms(st4,sr4)] do
                        for k3 in [1..NumberRealForms(st5,sr5)] do
                            for k4 in [1..NumberRealForms(st6,sr6)] do
                                H1:=RealFormById(st1,sr1,0);
                                H2:=RealFormById(st2,sr2,0);
                                H3:=RealFormById(st3,sr3,k1);
                                H4:=RealFormById(st4,sr4,k2);
                                H5:=RealFormById(st5,sr5,k3);
                                H6:=RealFormById(st6,sr6,k4);
                                rankRH:=RealRank(H1)+RealRank(H2)+RealRank(H3)+RealRank(H4)+RealRank(H5)+RealRank(H6);
                                dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2)+NonCompactDimension(H3)+NonCompactDimension(H4)+NonCompactDimension(H5)+NonCompactDimension(H6);
                                dimkH:=CompactDimension(H1)+CompactDimension(H2)+CompactDimension(H3)+CompactDimension(H4)+CompactDimension(H5)+CompactDimension(H6);   
                                if rankRG> rankRH  then
                                    if dimpG>= dimpH then
                                        if dimkG> dimkH then
                                            idH:=[rankRH,dimpH,st1,sr1,0,st2,sr2,0,st3,sr3,k1,st4,sr4,k2,st5,sr5,k3,st6,sr6,k4];
                                            Add(Sub,idH);
                                        fi;
                                    fi;
                                fi; 
                            od;
                        od;
                    od;
                od;
            fi;
            if Length(Separate(Separate(Separate(str)[2])[2]))=2 then
                s1:=Separate(str)[1];
                st1:=[s1[1]];
                sr1:=IntChar(s1[2])-48;
                s2:=Separate(Separate(str)[2])[1];
                st2:=[s2[1]];
                sr2:=IntChar(s2[2])-48;
                s3:=Separate(Separate(Separate(str)[2])[2])[1];
                st3:=[s3[1]];
                sr3:=IntChar(s3[2])-48;
                s4:=Separate(Separate(Separate(str)[2])[2])[2];
                st4:=[s4[1]];
                sr4:=IntChar(s4[2])-48;
                s5:=Separate(Separate(Separate(str)[2])[2])[2];
                st5:=[s5[4]];
                sr5:=IntChar(s5[5])-48;
                for k2 in [1..NumberRealForms(st4,sr4)] do
                    for k3 in [1..NumberRealForms(st5,sr5)] do
                        H1:=RealFormById(st1,sr1,0);
                        H2:=RealFormById(st2,sr2,0);
                        H3:=RealFormById(st3,sr3,0);
                        H4:=RealFormById(st4,sr4,k2);
                        H5:=RealFormById(st5,sr5,k3);
                        rankRH:=RealRank(H1)+RealRank(H2)+RealRank(H3)+RealRank(H4)+RealRank(H5);
                        dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2)+NonCompactDimension(H3)+NonCompactDimension(H4)+NonCompactDimension(H5);
                        dimkH:=CompactDimension(H1)+CompactDimension(H2)+CompactDimension(H3)+CompactDimension(H4)+CompactDimension(H5);   
                        if rankRG> rankRH  then
                            if dimpG>= dimpH then
                                if dimkG> dimkH then
                                    idH:=[rankRH,dimpH,st1,sr1,0,st2,sr2,0,st3,sr3,0,st4,sr4,k2,st5,sr5,k3];
                                    Add(Sub,idH);
                                fi;
                            fi;
                        fi; 
                    od;
                od;
            fi;
            if Length(Separate(Separate(Separate(Separate(str)[2])[2])[2]))=2 then
                s1:=Separate(str)[1];
                st1:=[s1[1]];
                sr1:=IntChar(s1[2])-48;
                s2:=Separate(Separate(str)[2])[1];
                st2:=[s2[1]];
                sr2:=IntChar(s2[2])-48;
                s3:=Separate(Separate(Separate(str)[2])[2])[1];
                st3:=[s3[1]];
                sr3:=IntChar(s3[2])-48;
                s4:=Separate(Separate(Separate(str)[2])[2])[2];
                st4:=[s4[1]];
                sr4:=IntChar(s4[2])-48;
                H1:=RealFormById(st1,sr1,0);
                H2:=RealFormById(st2,sr2,0);
                H3:=RealFormById(st3,sr3,0);
                H4:=RealFormById(st4,sr4,0);
                rankRH:=RealRank(H1)+RealRank(H2)+RealRank(H3)+RealRank(H4);
                dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2)+NonCompactDimension(H3)+NonCompactDimension(H4);
                dimkH:=CompactDimension(H1)+CompactDimension(H2)+CompactDimension(H3)+CompactDimension(H4);   
                if rankRG> rankRH  then
                    if dimpG>= dimpH then
                        if dimkG> dimkH then
                            idH:=[rankRH,dimpH,st1,sr1,0,st2,sr2,0,st3,sr3,0,st4,sr4,0];
                            Add(Sub,idH);
                        fi;
                    fi;
                fi; 
            fi;
        fi;
        if RemInt(i,10)=0 then
            Print("Counting: ");
            Print(String(Int(100*i/len)));
            Print("%");
            Print("\n");
        fi;        
    od;
    Print("Counting completed. \n");
    return Sub;
end);




###############################################################################
InstallGlobalFunction( PotentialSubalgebraPairs, function(arg)
local G,H,L,GG,type, rank, id, rankRG, dimpG,rankRH, dimpH,rankRL,dimpL,i,j,pairs,Hi,Lj ;
    type:=arg[1];
    rank:=arg[2];
    id:=arg[3];
    G:=RealFormById(type, rank, id);
    GG:=Concatenation(type,String(rank));
    rankRG:=RealRank(G);
    dimpG:=NonCompactDimension(G);
    H:=PotentialSubalgebras(type,rank,id);
    L:=H;
    pairs:=[];
    for i in [1..Length(H)] do
        for j in [1.. Length(L)] do
            Hi:=H[i];
            Lj:=L[j];
            rankRH:=Hi[1];
            dimpH:=Hi[2];
            rankRL:=Lj[1];
            dimpL:=Lj[2];
            if rankRG= rankRH+rankRL and dimpG=dimpH+dimpL and dimpH<= dimpL then
                Add(pairs,[Hi,Lj]);
            fi;
        od;
    od;
    return pairs;
end);

###############################################################################
InstallGlobalFunction( GetSymbolSimple, function(arg)
local type, rank, id, symb;
    type:=arg[1];
    rank:=arg[2];
    id:=arg[3];
    symb:=[];
    if type="A" then
        if rank=1 then
            if id=0 then
                symb:="sl(2,C)";
            elif id=1 then
                symb:="su(2)";
            elif id=2 then
                symb:="sl(2,R)";
            fi;
        elif rank=2 then
            if id=0 then
                symb:="sl(3,C)";
            elif id=1 then
                symb:="su(3)";
            elif id=2 then
                symb:="su(1,2)";
            elif id=3 then
                symb:="sl(3,R)";
            fi; 
        elif rank=3 then
            if id=0 then
                symb:="sl(4,C)";
            elif id=1 then
                symb:="su(4)";
            elif id=2 then
                symb:="su(1,3)";
            elif id=3 then
                symb:="su(2,2)";
            elif id=4 then
                symb:="su*(4)";
            elif id=5 then
                symb:="sl(4,R)";
            fi;     
        elif rank=4 then
            if id=0 then
                symb:="sl(5,C)";
            elif id=1 then
                symb:="su(5)";
            elif id=2 then
                symb:="su(1,4)";
            elif id=3 then
                symb:="su(2,3)";
            elif id=4 then
                symb:="sl(5,R)";
            fi;     
        elif rank=5 then
            if id=0 then
                symb:="sl(6,C)";
            elif id=1 then
                symb:="su(6)";
            elif id=2 then
                symb:="su(1,5)";
            elif id=3 then
                symb:="su(2,4)";
            elif id=4 then
                symb:="su(3,3)";
            elif id=5 then
                symb:="su*(6)";
            elif id=6 then
                symb:="sl(6,R)";
            fi;     
        elif rank=6 then
            if id=0 then
                symb:="sl(7,C)";
            elif id=1 then
                symb:="su(7)";
            elif id=2 then
                symb:="su(1,6)";
            elif id=3 then
                symb:="su(2,5)";
            elif id=4 then
                symb:="su(3,4)";
            elif id=5 then
                symb:="sl(7,R)";
            fi;
        elif rank=7 then
            if id=0 then
                symb:="sl(8,C)";
            elif id=1 then
                symb:="su(8)";
            elif id=2 then
                symb:="su(1,7)";
            elif id=3 then
                symb:="su(2,6)";
            elif id=4 then
                symb:="su(3,5)";
            elif id=5 then
                symb:="su(4,4)";
            elif id=6 then
                symb:="su*(8)";
            elif id=7 then
                symb:="sl(8,R)";
            fi;
        elif rank=8 then
            if id=0 then
                symb:="sl(9,C)";
            elif id=1 then
                symb:="su(9)";
            elif id=2 then
                symb:="su(1,8)";
            elif id=3 then
                symb:="su(2,7)";
            elif id=4 then
                symb:="su(3,6)";
            elif id=5 then
                symb:="su(4,5)";
            elif id=6 then
                symb:="sl(9,R)";
            fi;
        fi;
    elif type="B" then
        if rank=2 then
            if id=0 then
                symb:="so(5,C)";
            elif id=1 then
                symb:="so(5)";
            elif id=2 then
                symb:="so(2,3)";
            elif id=3 then
                symb:="so(4,1)";
            fi;
        elif rank=3 then
            if id=0 then
                symb:="so(7,C)";
            elif id=1 then
                symb:="so(7)";
            elif id=2 then
                symb:="so(2,5)";
            elif id=3 then
                symb:="so(4,3)";
            elif id=4 then
                symb:="so(6,1)";
            fi;
        elif rank=4 then
            if id=0 then
                symb:="so(9,C)";
            elif id=1 then
                symb:="so(9)";
            elif id=2 then
                symb:="so(2,7)";
            elif id=3 then
                symb:="so(4,5)";
            elif id=4 then
                symb:="so(6,3)";
            elif id=5 then
                symb:="so(8,1)";
            fi;
        elif rank=5 then
            if id=0 then
                symb:="so(11,C)";
            elif id=1 then
                symb:="so(11)";
            elif id=2 then
                symb:="so(2,9)";
            elif id=3 then
                symb:="so(4,7)";
            elif id=4 then
                symb:="so(6,5)";
            elif id=5 then
                symb:="so(8,3)";
            elif id=6 then
                symb:="so(10,1)";
            fi;
        elif rank=6 then
            if id=0 then
                symb:="so(13,C)";
            elif id=1 then
                symb:="so(13)";
            elif id=2 then
                symb:="so(2,11)";
            elif id=3 then
                symb:="so(4,9)";
            elif id=4 then
                symb:="so(6,7)";
            elif id=5 then
                symb:="so(8,5)";
            elif id=6 then
                symb:="so(10,3)";
            elif id=7 then
                symb:="so(12,1)";
            fi;
        elif rank=7 then
            if id=0 then
                symb:="so(15,C)";
            elif id=1 then
                symb:="so(15)";
            elif id=2 then
                symb:="so(2,13)";
            elif id=3 then
                symb:="so(4,11)";
            elif id=4 then
                symb:="so(6,9)";
            elif id=5 then
                symb:="so(8,7)";
            elif id=6 then
                symb:="so(10,5)";
            elif id=7 then
                symb:="so(12,3)";
            elif id=8 then
                symb:="so(14,1)";
            fi;
        elif rank=8 then
            if id=0 then
                symb:="so(17,C)";
            elif id=1 then
                symb:="so(17)";
            elif id=2 then
                symb:="so(2,15)";
            elif id=3 then
                symb:="so(4,13)";
            elif id=4 then
                symb:="so(6,11)";
            elif id=5 then
                symb:="so(8,9)";
            elif id=6 then
                symb:="so(10,7)";
            elif id=7 then
                symb:="so(12,5)";
            elif id=8 then
                symb:="so(14,3)";
            elif id=9 then
                symb:="so(16,1)";
            fi;
        fi;
    elif type="C" then
        if rank=3 then
            if id=0 then
                symb:="sp(3,C)";
            elif id=1 then
                symb:="sp(3)";
            elif id=2 then
                symb:="sp(1,2)";
            elif id=3 then
                symb:="sp(3,R)";
            fi;
        elif rank=4 then
            if id=0 then
                symb:="sp(4,C)";
            elif id=1 then
                symb:="sp(4)";
            elif id=2 then
                symb:="sp(1,3)";
            elif id=3 then
                symb:="sp(2,2)";
            elif id=4 then
                symb:="sp(4,R)";
            fi;
        elif rank=5 then
            if id=0 then
                symb:="sp(5,C)";
            elif id=1 then
                symb:="sp(5)";
            elif id=2 then
                symb:="sp(1,4)";
            elif id=3 then
                symb:="sp(2,3)";
            elif id=4 then
                symb:="sp(5,R)";
            fi;
        elif rank=6 then
            if id=0 then
                symb:="sp(6,C)";
            elif id=1 then
                symb:="sp(6)";
            elif id=2 then
                symb:="sp(1,5)";
            elif id=3 then
                symb:="sp(2,4)";
            elif id=4 then
                symb:="sp(3,3)";
            elif id=5 then
                symb:="sp(6,R)";
            fi;
        elif rank=7 then
            if id=0 then
                symb:="sp(7,C)";
            elif id=1 then
                symb:="sp(7)";
            elif id=2 then
                symb:="sp(1,6)";
            elif id=3 then
                symb:="sp(2,5)";
            elif id=4 then
                symb:="sp(3,4)";
            elif id=5 then
                symb:="sp(7,R)";
            fi;
        elif rank=8 then
            if id=0 then
                symb:="sp(8,C)";
            elif id=1 then
                symb:="sp(8)";
            elif id=2 then
                symb:="sp(1,7)";
            elif id=3 then
                symb:="sp(2,6)";
            elif id=4 then
                symb:="sp(3,5)";
            elif id=5 then
                symb:="sp(4,4)";
            elif id=6 then
                symb:="sp(8,R)";
            fi;
        fi;
    elif type="D" then
        if rank=4 then
            if id=0 then
                symb:="so(8,C)";
            elif id=1 then
                symb:="so(8)";
            elif id=2 then
                symb:="so(2,6)";
            elif id=3 then
                symb:="so(4,4)";
            elif id=4 then
                symb:="so(3,5)";
            elif id=5 then
                symb:="so(1,7)";
            fi;
        elif rank=5 then
            if id=0 then
                symb:="so(10,C)";
            elif id=1 then
                symb:="so(10)";
            elif id=2 then
                symb:="so(2,8)";
            elif id=3 then
                symb:="so(4,6)";
            elif id=4 then
                symb:="so*(10)";
            elif id=5 then
                symb:="so(9,1)";
            elif id=6 then
                symb:="so(3,7)";
            elif id=7 then
                symb:="so(5,5)";
            fi;
        elif rank=6 then
            if id=0 then
                symb:="so(12,C)";
            elif id=1 then
                symb:="so(12)";
            elif id=2 then
                symb:="so(2,10)";
            elif id=3 then
                symb:="so(4,8)";
            elif id=4 then
                symb:="so(6,6)";
            elif id=5 then
                symb:="so*(12)";
            elif id=6 then
                symb:="so(11,1)";
            elif id=7 then
                symb:="so(3,9)";
            elif id=8 then
                symb:="so(5,7)";
            fi;
        elif rank=7 then
            if id=0 then
                symb:="so(14,C)";
            elif id=1 then
                symb:="so(14)";
            elif id=2 then
                symb:="so(2,12)";
            elif id=3 then
                symb:="so(4,10)";
            elif id=4 then
                symb:="so(6,8)";
            elif id=5 then
                symb:="so*(14)";
            elif id=6 then
                symb:="so(13,1)";
            elif id=7 then
                symb:="so(3,11)";
            elif id=8 then
                symb:="so(5,9)";
            elif id=9 then
                symb:="so(7,7)";
            fi;
        elif rank=8 then
            if id=0 then
                symb:="so(16,C)";
            elif id=1 then
                symb:="so(16)";
            elif id=2 then
                symb:="so(2,14)";
            elif id=3 then
                symb:="so(4,12)";
            elif id=4 then
                symb:="so(6,10)";
            elif id=5 then
                symb:="so(8,8)";
            elif id=6 then
                symb:="so*(16)";
            elif id=7 then
                symb:="so(15,1)";
            elif id=8 then
                symb:="so(3,13)";
            elif id=9 then
                symb:="so(5,11)";
            elif id=10 then
                symb:="so(7,9)";
            fi;
        fi;
    elif type="E" then
        if rank=6 then
            if id=0 then
                symb:="E(6,C)";
            elif id=1 then
                symb:="E6(-78)";
            elif id=2 then
                symb:="E6(6)";
            elif id=3 then
                symb:="E6(2)";
            elif id=4 then
                symb:="E6(-14)";
            elif id=5 then
                symb:="E6(-26)";
            fi;
        elif rank=7 then
            if id=0 then
                symb:="E(7,C)";
            elif id=1 then
                symb:="E7(-133)";
            elif id=2 then
                symb:="E7(7)";
            elif id=3 then
                symb:="E7(-25)";
            elif id=4 then
                symb:="E7(-5)";
            fi;
        elif rank=8 then
            if id=0 then
                symb:="E(8,C)";
            elif id=1 then
                symb:="E8(-248)";
            elif id=2 then
                symb:="E8(8)";
            elif id=3 then
                symb:="E8(-24)";
            fi;
        fi;
    elif type="F" then
        if rank=4 then
            if id=0 then
                symb:="F(4,C)";
            elif id=1 then
                symb:="F4(-52)";
            elif id=2 then
                symb:="F4(4)";
            elif id=3 then
                symb:="F4(-20)";
            fi;
        fi;
    elif type="G" then
        if rank=2 then
            if id=0 then
                symb:="G(2,C)";
            elif id=1 then
                symb:="G2(-14)";
            elif id=2 then
                symb:="G2(2)";
            fi;
        fi;
    fi;
    return symb;
end);

###############################################################################
InstallGlobalFunction( GetSymbolSemisimple, function(arg)
local t1, symb;
    t1:=arg[1];
    symb:=GetSymbolSimple(t1[3],t1[4],t1[5]);
    if Length(t1)>5 then
        symb:=Concatenation(symb,"+",GetSymbolSimple(t1[6],t1[7],t1[8]));
    fi;
    if Length(t1)>8 then
        symb:=Concatenation(symb,"+",GetSymbolSimple(t1[9],t1[10],t1[11]));
    fi;
    if Length(t1)>11 then
        symb:=Concatenation(symb,"+",GetSymbolSimple(t1[12],t1[13],t1[14]));
    fi;
    if Length(t1)>14 then
        symb:=Concatenation(symb,"+",GetSymbolSimple(t1[15],t1[16],t1[17]));
    fi;
    if Length(t1)>17 then
        symb:=Concatenation(symb,"+",GetSymbolSimple(t1[18],t1[19],t1[20]));
    fi;
    if Length(t1)>20 then
        symb:=Concatenation(symb,"+",GetSymbolSimple(t1[21],t1[22],t1[23]));
    fi;
    if Length(t1)>23 then
        symb:=Concatenation(symb,"+",GetSymbolSimple(t1[24],t1[25],t1[26]));
    fi;
    return symb;
end);

###############################################################################
InstallGlobalFunction( GetSymbolPairs, function(arg)
local input, i,n,output, pairs,t1,t2;
    input:=arg[1];
    n:=Length(input);
    if n>0 then
        for i in [1..n] do
            pairs:=input[i];
            t1:=pairs[1];
            t2:=pairs[2];
            Print("h=");
            Print(GetSymbolSemisimple(t1));
            Print(", ");
            Print("l=");
            Print(GetSymbolSemisimple(t2));
            Print(" \n");
        od;
    fi;
end);

###############################################################################
InstallGlobalFunction( RealFormByTuple, function(arg)
local t, G,H;
    t:=arg[1];
    G:=RealFormById(t[3],t[4],t[5]);
    if Length(t)>5 then
        H:=RealFormById(t[6],t[7],t[8]);
        G:=DirectSumOfAlgebras(G,H);
    fi;
    if Length(t)>8 then
        H:=RealFormById(t[9],t[10],t[11]);
        G:=DirectSumOfAlgebras(G,H);
    fi;
    if Length(t)>11 then
        H:=RealFormById(t[12],t[13],t[14]);
        G:=DirectSumOfAlgebras(G,H);
    fi;
    if Length(t)>14 then
        H:=RealFormById(t[15],t[16],t[17]);
        G:=DirectSumOfAlgebras(G,H);
    fi;
    if Length(t)>17 then
        H:=RealFormById(t[18],t[19],t[20]);
        G:=DirectSumOfAlgebras(G,H);
    fi;
    if Length(t)>20 then
        H:=RealFormById(t[21],t[22],t[23]);
        G:=DirectSumOfAlgebras(G,H);
    fi;
    if Length(t)>23 then
        H:=RealFormById(t[24],t[25],t[26]);
        G:=DirectSumOfAlgebras(G,H);
    fi;
    return G;
end);

###############################################################################
InstallGlobalFunction( CheckTuple, function(arg)
local t, log, G;
    t:=arg[1];
    log:=false;
    G:=RealFormByTuple(t);
    if RealRank(G)=t[1] and NonCompactDimension(G)=t[2] then
        log:=true;
    fi;
    return log;
end);


###############################################################################
InstallGlobalFunction( CheckRankConditions, function(arg)
local typeg, rankg, idg,maxs,G,H,l0,l1,l2,l3, index;
    typeg:=arg[1];
    rankg:=arg[2];
    idg:=arg[3];
    G:=RealFormById(typeg,rankg,idg);
    maxs:=MaximalReductiveSubalgebras( typeg, rankg, idg );
    Print("g=",NameRealForm(G)); 
    Print(" | real rank(g)=",RealRank(G)) ;  
    Print(" | a-hyp rank(g)=",AHypRank(G),"\n") ;
    Print("----------------------------\n");
#   Print("L0-Calabi–Markus phenomenon real rank(g)=real rank(h)\n");  
#   Print("L1-ahyp rank(g)=ahyp rank(h)\n");  
#   Print("L2-ahyp rank(g)>real rank(h)\n");  
#   Print("L3-none of the above conditions is met\n");  
#   Print("----------------------------\n");
    index:=1;
    for H in maxs.subalgs do
        Print("#",index,": h=",NameRealForm(H)," | real rank(h)=",RealRank(H)," | ahyp rank(h)=",AHypRank(H),"\n");
        l0:=RealRank(G)=RealRank(H);
        Print(" | L0-",l0);
        l1:=AHypRank(G)=AHypRank(H);
        Print(" | L1-",l1);
        l2:=AHypRank(G)>RealRank(H);
        Print(" | L2-",l2);
        l3:=not(l0 or l1 or l2);
        Print(" | L3-",l3);
        Print("\n");
        Print("----\n");
        index:=index+1;
    od;
end);

###############################################################################
InstallGlobalFunction( TechCoeff, function(arg)
local jgc,jhc,out, out2,k,z2,l,w2,v2,j;
    jgc:=arg[1];
    jhc:=arg[2];
    out:=[];
    for j in [1..Length(Basis(jhc))] do
        out2:=ListWithIdenticalEntries( Length(Basis(jgc)), 0 );;
        v2:=String(Basis(jhc)[j]);
        w2:=SplitString(v2,"+");
        for k in [1..Length(w2)] do
            z2:=w2[k];
            z2:=SplitString(z2,"*");
            for l in [1..Length(Basis(jgc))] do
                if Length(z2)=1 and z2[1]=String(Basis(jgc)[l]) then
                    out2[l]:=1;
                elif Length(z2)=2 and z2[2]=String(Basis(jgc)[l]) then
                    out2[l]:=EvalString(z2[1]);
                fi;
            od;
        od;
        Add(out,out2);
    od;
    return out;
end);


###############################################################################
InstallGlobalFunction( CheckProperSL2RAction, function(arg)
local idG,G,Hmax, idH,nHmax, idhc,j, idhc2,j2,GC,orbs,jg,B,W,ki,nb,nb2,s,sub,g,m,lo,i,subs,res,el,k,HC,jhc,jgc,hs,cf,orb,k2,res2,w,li;
    idG:=[arg[1],arg[2],arg[3]];
    Hmax:=MaximalReductiveSubalgebras(idG[1],idG[2],idG[3]).subalgs[arg[4]];
    nHmax:=NameRealForm(Hmax);
    nHmax:=SplitString(nHmax,"+");
    idhc:=[];
    for j in [1..Length(nHmax)] do
        Add(idhc,GetIdFromName(nHmax[j]));
    od;
    idhc2:="";
    for j in [1..Length(idhc)] do
        if j>1 then
            idhc2:=Concatenation(idhc2," ");
        fi;
        j2:=idhc[j];
        idhc2:=Concatenation(idhc2,j2[1],String(j2[2]));
    od;
    GC:= SimpleLieAlgebra( idG[1],idG[2], Rationals );
    orbs:= NilpotentOrbits(GC);;
    jg:=CartanSubalgebra(GC);;
    B:=Basis(GC);;
    W:= WeylGroup( RootSystem( GC ) );;
    lo:=Length(orbs);;
    nb:=Length(B);;
    nb2:=nb-idG[2]+1;;
    s:= LieAlgebraAndSubalgebras( Concatenation(idG[1],String(idG[2])) );;
    sub:= s.subalgs;;
    g:= InclusionsGraph( Concatenation(idG[1],String(idG[2])) );;
    m:= Filtered( g, x -> x[1]=0 );; i:= List( m, x -> x[2] );
    subs:=List( sub{i}, SemiSimpleType );
    res2:=[];
    for k  in [1..Length(i)] do
        if SemiSimpleType(sub[i[k]])= idhc2 then
            Add(res2,i[k]);
        fi;
    od;
    for k2 in [1..Length(res2)] do
        HC:=sub[res2[k2]];
        jhc:=CartanSubalgebra(HC);
        jgc:=CartanSubalgebra(GC);
        for i in [1..lo] do
            hs:=SL2Triple( orbs[i] )[2];
            cf:=Coefficients( B,hs){[nb2..nb]};
            orb:= WeylOrbitIterator( W, cf );
            k:=0;
            while not IsDoneIterator( orb ) do
                res:=TechCoeff(jgc,jhc);
                w:= NextIterator( orb ); Add(res,List(w));
                li:=RankMat(res);
                if li <Length(res) then
                    k:=k+1;
                fi; 
            od;
            if k =0 then
                Print("proper\n");
                return;
            fi;
        od;
    od;
    Print("not proper\n");
end);

#E  ckforms.gi  . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here