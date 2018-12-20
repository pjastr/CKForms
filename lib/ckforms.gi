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
InstallGlobalFunction( CompactDimension, function(arg)
local G, K, dimK;
    G:=arg[1];
    K:=CartanDecomposition( G ).K;
    dimK:= Dimension(K);
    return dimK;
end);

###############################################################################
InstallGlobalFunction( Simple, function()
local ListSimple,i,k,j ;
    ListSimple:=[];
    #A
    Print("A \n");
    for k in [1..8] do
        for i in [1..NumberRealForms("A",k)] do
            Append(ListSimple,[["A",k,i]]);
        od;
    od;
    #B
    Print("B \n");
    for k in [2..8] do
        for i in [1..NumberRealForms("B",k)] do
            Append(ListSimple,[["B",k,i]]);
        od;
    od;
    #C
    Print("C \n");
    for k in [3..8] do
        for i in [1..NumberRealForms("C",k)] do
            Append(ListSimple,[["C",k,i]]);
        od;
    od;
    #D
    Print("D \n");
    for k in [4..8] do
        for i in [1..NumberRealForms("D",k)] do
            Append(ListSimple,[["D",k,i]]);
        od;
    od;
    #E
    Print("E \n");
    for k in [6..8] do
        for i in [1..NumberRealForms("E",k)] do
            Append(ListSimple,[["E",k,i]]);
        od;
    od;
    #F
    Print("F \n");
    for i in [1..NumberRealForms("F",4)] do
        Append(ListSimple,[["F",4,i]]);
    od;
    #G
    Print("G \n");
    for i in [1..NumberRealForms("G",2)] do
        Append(ListSimple,[["G",2,i]]);
    od;
    return ListSimple;
end);

###############################################################################
InstallGlobalFunction( Simple2, function()
local ListSR,k,A,i,r1,r2, ListSR2;
    A:=Simple();
    ListSR:=Tuples(A,2);
    Print("realification \n");
    #A
    for k in [1..4] do
        Append(ListSR,[[["A",k,0],["A",k,0]]]);
    od;
    #B
    for k in [2..4] do
        Append(ListSR,[[["B",k,0],["B",k,0]]]);
    od;
    #C
    for k in [3..4] do
        Append(ListSR,[[["C",k,0],["C",k,0]]]);
    od;
    #D
    Append(ListSR,[[["D",4,0],["D",4,0]]]);
    #E - none realification with rank<=8
    #F
    Append(ListSR,[[["F",4,0],["F",4,0]]]);
    #G
    Append(ListSR,[[["G",2,0],["G",2,0]]]);
    ListSR2:=[];
    for i in [1..Length(ListSR)] do
        r1:=ListSR[i][1][2];
        r2:=ListSR[i][2][2];
        if r1+r2<=8 then
            Append(ListSR2,[ListSR[i]]);
        fi;
    od;
    return ListSR2;
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
    ListSub:=Immutable(ListSub);
    return ListSub;
end);



###############################################################################
InstallGlobalFunction( GetAllForms, function(arg)
local G, AF, len,i,j, forms,t1,r1,t2,r2;
    G:=arg[1];
    len:=Length(G);
    AF:=[];
    t1:=[G[1]];
    r1:=IntChar(G[2])-48;
    if len>2 then
        t2:=[G[4]];
        r2:=IntChar(G[5])-48;
    fi;
    if len=2 then
        for i in [1..NumberRealForms(t1,r1)] do
            Append(AF,[[[t1,r1,i]]]);
        od;
    elif len=5 then
        for i in [1..NumberRealForms(t1,r1)] do
            for j in [1..NumberRealForms(t2,r2)] do
                Append(AF,[[[t1,r1,i],[t2,r2,j]]]);
            od;
        od;
        if t1=t2 and r1=r2 then
            Append(AF,[[[t1,r1,0],[t1,r1,0]]]);
        fi;
    fi;
    return AF;
end);

###############################################################################
InstallGlobalFunction( PotentialSubalgebras, function(arg)
local G,H,GG,type, rank, id, str,len,t1,r1,t2,r2,t3,r3,i,j,k,k1, ListSub,Sub, rankRG, rankRH, dimpG, dimkG;
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
        if Length(str)>2 then
            t2:=[str[4]];
            r2:=IntChar(str[5])-48;
        fi;
        if Length(str)>5 then
            t3:=[str[7]];
            r3:=IntChar(str[8])-48;
        fi;
        Print("\n step1");
        if Length(str) =2 then
            for k in [1..NumberRealForms(t1,r1)] do
                H:=RealFormById(t1,r1,k);
                if rankRG> RealRank(H)  then
                    if dimpG>= NonCompactDimension(H) then
                        if dimkG> CompactDimension(H) then
                            Add(Sub,H);
                        fi;
                    fi;
                fi;
            od;
        fi;
        Print("step2");
        if Length(str) =5 then
            for k in [1..NumberRealForms(t1,r1)] do
                for j in [1..NumberRealForms(t2,r2)] do
                    H:=DirectSumOfAlgebras(RealFormById(t1,r1,k),RealFormById(t2,r2,j));
                    if rankRG> RealRank(H)  then
                        if dimpG>= NonCompactDimension(H) then
                            if dimkG> CompactDimension(H) then
                                Add(Sub,H);
                            fi;
                        fi;
                    fi;
                od;
            od;
            Print("step2a");
            if t1=t2 and r1=r2 then
                H:=RealFormById(t1,r1,0);
                if rankRG> RealRank(H)  then
                    if dimpG>= NonCompactDimension(H) then
                        if dimkG> CompactDimension(H) then
                            Add(Sub,H);
                        fi;
                    fi;
                fi;
            fi;
        fi;
        Print("step3");
        if Length(str) =8 then
            for k in [1..NumberRealForms(t1,r1)] do
                for j in [1..NumberRealForms(t2,r2)] do
                    for k1 in [1..NumberRealForms(t3,r3)] do
                        H:=DirectSumOfAlgebras(DirectSumOfAlgebras(RealFormById(t1,r1,k),RealFormById(t2,r2,j)),RealFormById(t3,r3,k1));
                        if rankRG> RealRank(H)  then
                            if dimpG>= NonCompactDimension(H) then
                                if dimkG> CompactDimension(H) then
                                    Add(Sub,H);
                                fi;
                            fi;
                        fi;
                    od;
                od;
            od;
            # if t1=t2 and r1=r2 then
            #     H:=RealFormById(t1,r1,0);
            #     if rankRG> RealRank(H) and dimpG> NonCompactDimension(H) and dimkG> CompactDimension(H) then
            #         Add(Sub,H);
            #     fi;
            # fi;
        fi;
        Print(String(i));
    od;
    return Sub;
end);

###############################################################################
InstallGlobalFunction( NicerPotentialSubalgebras, function(arg)
local G,H1,H2,H3,H4,H5,H6,H7,H8,GG,type, rank, id,idH, str,len,t1,r1,t2,r2,t3,r3,t4,r4,t5,r5,t6,r6,t7,r7,t8,r8,i,j,k,k1,k2,k3,k4,k5,k6,ListSub,Sub, rankRG, rankRH, dimpG, dimpH, dimkG, dimkH;
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
            if t1=t2 and r1=r2 then
                H1:=RealFormById(t1,r1,0);
                rankRH:=RealRank(H1);
                dimpH:=NonCompactDimension(H1);
                dimkH:=CompactDimension(H1);
                if rankRG> rankRH  then
                    if dimpG>= dimpH then
                        if dimkG> dimkH then
                            idH:=[rankRH,dimpH,t1,r1,0,t1,r1,0];
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
            if t1=t2 and r1=r2 then
                for k1 in [1..NumberRealForms(t3,r3)] do
                    H1:=RealFormById(t1,r1,0);
                    H2:=RealFormById(t3,r3,k1);
                    rankRH:=RealRank(H1)+RealRank(H2);
                    dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2);
                    dimkH:=CompactDimension(H1)+CompactDimension(H2);  
                    if rankRG> rankRH  then
                        if dimpG>= dimpH then
                            if dimkG> dimkH then
                                idH:=[rankRH,dimpH,t1,r1,0,t1,r1,0,t3,r3,k1];
                                Add(Sub,idH);
                            fi;
                        fi;
                    fi;
                od;
            fi;
            if t2=t3 and r2=r3 then
                for k in [1..NumberRealForms(t1,r1)] do
                    H1:=RealFormById(t1,r1,k);
                    H2:=RealFormById(t2,r2,0);
                    rankRH:=RealRank(H1)+RealRank(H2);
                    dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2);
                    dimkH:=CompactDimension(H1)+CompactDimension(H2);  
                    if rankRG> rankRH  then
                        if dimpG>= dimpH then
                            if dimkG> dimkH then
                                idH:=[rankRH,dimpH,t1,r1,k,t2,r2,0,t2,r2,0];
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
            if t1=t2 and r1=r2 then
                for k1 in [1..NumberRealForms(t3,r3)] do
                    for k2 in [1..NumberRealForms(t4,r4)] do
                        H1:=RealFormById(t1,r1,0);
                        H2:=RealFormById(t3,r3,k1);
                        H3:=RealFormById(t4,r4,k2);
                        rankRH:=RealRank(H1)+RealRank(H2)+RealRank(H3);
                        dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2)+NonCompactDimension(H3);
                        dimkH:=CompactDimension(H1)+CompactDimension(H2)+CompactDimension(H3);  
                        if rankRG> rankRH  then
                            if dimpG>= dimpH then
                                if dimkG> dimkH then
                                    idH:=[rankRH,dimpH,t1,r1,0,t1,r1,0,t3,r3,k1,t4,r4,k2];
                                    Add(Sub,idH);
                                fi;
                            fi;
                        fi;
                    od;
                od;
            fi;
            if t2=t3 and r2=r3 then
                for k in [1..NumberRealForms(t1,r1)] do
                    for k2 in [1..NumberRealForms(t4,r4)] do
                        H1:=RealFormById(t1,r1,k);
                        H2:=RealFormById(t2,r2,0);
                        H3:=RealFormById(t4,r4,k2);
                        rankRH:=RealRank(H1)+RealRank(H2)+RealRank(H3);
                        dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2)+NonCompactDimension(H3);
                        dimkH:=CompactDimension(H1)+CompactDimension(H2)+CompactDimension(H3);  
                        if rankRG> rankRH  then
                            if dimpG>= dimpH then
                                if dimkG> dimkH then
                                    idH:=[rankRH,dimpH,t1,r1,k,t2,r2,0,t2,r2,0,t4,r4,k2];
                                    Add(Sub,idH);
                                fi;
                            fi;
                        fi;
                    od;
                od;
            fi;
            if t3=t4 and r3=r4 then
                for k in [1..NumberRealForms(t1,r1)] do
                    for j in [1..NumberRealForms(t2,r2)] do
                        H1:=RealFormById(t1,r1,k);
                        H2:=RealFormById(t2,r2,j);
                        H3:=RealFormById(t3,r3,0);
                        rankRH:=RealRank(H1)+RealRank(H2)+RealRank(H3);
                        dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2)+NonCompactDimension(H3);
                        dimkH:=CompactDimension(H1)+CompactDimension(H2)+CompactDimension(H3);  
                        if rankRG> rankRH  then
                            if dimpG>= dimpH then
                                if dimkG> dimkH then
                                    idH:=[rankRH,dimpH,t1,r1,k,t2,r2,j,t3,r3,0,t3,r3,0];
                                    Add(Sub,idH);
                                fi;
                            fi;
                        fi;
                    od;
                od;
            fi;
            if t1=t2 and r1=r2 and t3=t4 and r3=r4 then            
                H1:=RealFormById(t1,r1,0);
                H2:=RealFormById(t3,r3,0);
                rankRH:=RealRank(H1)+RealRank(H2);
                dimpH:=NonCompactDimension(H1)+NonCompactDimension(H2);
                dimkH:=CompactDimension(H1)+CompactDimension(H2);  
                if rankRG> rankRH  then
                    if dimpG>= dimpH then
                        if dimkG> dimkH then
                            idH:=[rankRH,dimpH,t1,r1,0,t2,r2,0,t3,r3,0,t4,r4,0];
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
        fi;
        Print(String(i));
        Print("\n");
    od;
    return Sub;
end);




###############################################################################
InstallGlobalFunction( NicerPotentialTriples, function(arg)
local G,H,L,GG,type, rank, id, rankRG, dimpG,rankRH, dimpH,rankRL,dimpL,i,j,Triples,Hi,Lj ;
    type:=arg[1];
    rank:=arg[2];
    id:=arg[3];
    G:=RealFormById(type, rank, id);
    GG:=Concatenation(type,String(rank));
    rankRG:=RealRank(G);
    dimpG:=NonCompactDimension(G);
    H:=NicerPotentialSubalgebras(type,rank,id);
    L:=H;
    Triples:=[];
    for i in [1..Length(H)] do
        for j in [1.. Length(L)] do
            Hi:=H[i];
            Lj:=L[j];
            rankRH:=Hi[1];
            dimpH:=Hi[2];
            rankRL:=Lj[1];
            dimpL:=Lj[2];
            if rankRG= rankRH+rankRL and dimpG=dimpH+dimpL and dimpH<= dimpL then
                Add(Triples,[Hi,Lj]);
            fi;
        od;
    od;
    return Triples;
end);

#E  ckforms.gi  . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here