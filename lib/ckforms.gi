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
    ListSub:=Immutable(ListSub);
    return ListSub;
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
local G,H1,H2,H3,H4,H5,H6,H7,H8,GG,type, rank, id,idH, str,len,t1,r1,t2,r2,t3,r3,t4,r4,t5,r5,t6,r6,t7,r7,t8,r8,i,j,k,k1,k2,k3,k4,k5,k6,ListSub,Sub, rankRG, rankRH, dimpG, dimpH, dimkG, dimkH, s1, st1, sr1,s2,st2,sr2,s3,st3,sr3,s4,st4,sr4;
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
                            idH:=[rankRH,dimpH,st1,sr1,0,st1,sr1,0];
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
                                idH:=[rankRH,dimpH,st1,sr1,0,st1,sr1,0,st2,sr2,j];
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
                                    idH:=[rankRH,dimpH,st1,sr1,0,st1,sr1,0,st2,sr2,j,st3,sr3,k1];
                                    Add(Sub,idH);
                                fi;
                            fi;
                        fi; 
                    od;
                od;
            fi;
            if Length(Separate(Separate(str)[2]))>2 then
                s1:=Separate(str)[1];
                st1:=[s1[1]];
                sr1:=IntChar(s1[2])-48;
                s2:=Separate(s1)[2];
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
                            idH:=[rankRH,dimpH,st1,sr1,0,st1,sr1,0,st2,sr2,0,st2,sr2,0];
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
                                        idH:=[rankRH,dimpH,st1,sr1,0,st1,sr1,0,st2,sr2,j,st3,sr3,k1,st4,sr4,k2];
                                        Add(Sub,idH);
                                    fi;
                                fi;
                            fi; 
                        od;
                    od;
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

###############################################################################
InstallGlobalFunction( GetSymbolRealForm, function(arg)
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
                symb:="sl(2,H)=su*(4)";
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
                symb:="sl(3,H)=su*(6)";
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
                symb:="sl(4,H)=su*(8)";
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
                symb:="sp(3,C)/sp(6,C)";
            elif id=1 then
                symb:="sp(3)";
            elif id=2 then
                symb:="sp(1,2)";
            elif id=3 then
                symb:="sp(3,R)/sp(6,R)";
            fi;
        elif rank=4 then
            if id=0 then
                symb:="sp(4,C)/sp(8,C)";
            elif id=1 then
                symb:="sp(4)";
            elif id=2 then
                symb:="sp(1,3)";
            elif id=3 then
                symb:="sp(2,2)";
            elif id=4 then
                symb:="sp(4,R)/sp(8,R)";
            fi;
        elif rank=5 then
            if id=0 then
                symb:="sp(5,C)/sp(10,C)";
            elif id=1 then
                symb:="sp(5)";
            elif id=2 then
                symb:="sp(1,4)";
            elif id=3 then
                symb:="sp(2,3)";
            elif id=4 then
                symb:="sp(5,R)/sp(10,R)";
            fi;
        elif rank=6 then
            if id=0 then
                symb:="sp(6,C)/sp(12,C)";
            elif id=1 then
                symb:="sp(6)";
            elif id=2 then
                symb:="sp(1,5)";
            elif id=3 then
                symb:="sp(2,4)";
            elif id=4 then
                symb:="sp(3,3)";
            elif id=5 then
                symb:="sp(6,R)/sp(12,R)";
            fi;
        elif rank=7 then
            if id=0 then
                symb:="sp(7,C)/sp(14,C)";
            elif id=1 then
                symb:="sp(7)";
            elif id=2 then
                symb:="sp(1,6)";
            elif id=3 then
                symb:="sp(2,5)";
            elif id=4 then
                symb:="sp(3,4)";
            elif id=5 then
                symb:="sp(7,R)/sp(14,R)";
            fi;
        elif rank=8 then
            if id=0 then
                symb:="sp(8,C)/sp(16,C)";
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
                symb:="sp(8,R)/sp(16,R)";
            fi;
        fi;
    elif type="D" then
        if rank=4 then
            if id=0 then
                symb:="so(8,C)";
            elif id=1 then
                symb:="so(8)";
            elif id=2 then
                symb:="so*(8)=so(2,6)";
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
                symb:="EI=E6(6)";
            elif id=3 then
                symb:="EII=E6(2)";
            elif id=4 then
                symb:="EIII=E6(-14)";
            elif id=5 then
                symb:="EIV=E6(-26)";
            fi;
        elif rank=7 then
            if id=0 then
                symb:="E(7,C)";
            elif id=1 then
                symb:="E7(-133)";
            elif id=2 then
                symb:="EV=E7(7)";
            elif id=3 then
                symb:="EVII=E7(-25)";
            elif id=4 then
                symb:="EVI=E7(-5)";
            fi;
        elif rank=8 then
            if id=0 then
                symb:="E(8,C)";
            elif id=1 then
                symb:="E8(-248)";
            elif id=2 then
                symb:="EVIII=E8(8)";
            elif id=3 then
                symb:="EIX=E8(-24)";
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
InstallGlobalFunction( GetSymbolTriples, function(arg)
local input, i,n,output, triple,t1,t2;
    input:=arg[1];
    n:=Length(input);
    if n>0 then
        for i in [1..n] do
            triple:=input[i];
            t1:=triple[1];
            t2:=triple[2];
            Print("h=");
            Print(GetSymbolRealForm(t1[3],t1[4],t1[5]));
            if (Length(t1))>5 then
                Print("+");
                Print(GetSymbolRealForm(t1[6],t1[7],t1[8]));
            fi;
            if (Length(t1))>8 then
                Print("+");
                Print(GetSymbolRealForm(t1[9],t1[10],t1[11]));
            fi;
            if (Length(t1))>11 then
                Print("+");
                Print(GetSymbolRealForm(t1[12],t1[13],t1[14]));
            fi;
            if (Length(t1))>14 then
                Print("+");
                Print(GetSymbolRealForm(t1[15],t1[16],t1[17]));
            fi;
            if (Length(t1))>17 then
                Print("+");
                Print(GetSymbolRealForm(t1[18],t1[19],t1[20]));
            fi;
            if (Length(t1))>20 then
                Print("+");
                Print(GetSymbolRealForm(t1[21],t1[22],t1[23]));
            fi;
            if (Length(t1))>23 then
                Print("+");
                Print(GetSymbolRealForm(t1[24],t1[25],t1[26]));
            fi;
            Print(", ");
            Print("l=");
            Print(GetSymbolRealForm(t2[3],t2[4],t2[5]));
            if (Length(t2))>5 then
                Print("+");
                Print(GetSymbolRealForm(t2[6],t2[7],t2[8]));
            fi;
            if (Length(t2))>8 then
                Print("+");
                Print(GetSymbolRealForm(t2[9],t2[10],t2[11]));
            fi;
            if (Length(t2))>11 then
                Print("+");
                Print(GetSymbolRealForm(t2[12],t2[13],t2[14]));
            fi;
            if (Length(t2))>14 then
                Print("+");
                Print(GetSymbolRealForm(t2[15],t2[16],t2[17]));
            fi;
            if (Length(t2))>17 then
                Print("+");
                Print(GetSymbolRealForm(t2[18],t2[19],t2[20]));
            fi;
            if (Length(t2))>20 then
                Print("+");
                Print(GetSymbolRealForm(t2[21],t2[22],t2[23]));
            fi;
            if (Length(t2))>23 then
                Print("+");
                Print(GetSymbolRealForm(t2[24],t2[25],t2[26]));
            fi;
            Print(" \n");
        od;
    fi;
end);

#E  ckforms.gi  . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here