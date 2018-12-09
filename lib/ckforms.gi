#############################################################################
##
#W  ckforms.gi                  CKForms Package                  
##
##  Installation file for functions of the CKForms package.
##


###############################################################################
InstallGlobalFunction( NonCompactDimension, function(arg)
local G, dimG, K, dimK, dG;
    G:=arg[1];
    dimG:=Dimension(G);
    K:=CartanDecomposition( G ).K;
    dimK:= Dimension(K);
    dG:=dimG-dimK;
    return dG;
end);

#E  ckforms.gi  . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here