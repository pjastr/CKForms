#############################################################################
##
##
#W  ckforms.gd                   CKForms Package                  
##
##  Declaration file for functions of the CKForms package.
##



#############################################################################
##
#F  NonCompactDimension(g)
##
##  The input is a real Lie algebra g. This function returns non-compact dimension
##  of g (dimension of non-compact part from the Cartan decomposition for g).
##
DeclareGlobalFunction( "NonCompactDimension" );

#############################################################################
##
#F  RealRank(g)
##
##  The input is a real Lie algebra g (as an object, not as a tuple). The output
##  is the real rank of g (the dimension of the Cartan subalgebra of g).
##
DeclareGlobalFunction( "RealRank" );

#############################################################################
##
#F  CompactDimension(g)
##
##  The input is a real Lie algebra g. This function returns the compact dimension
##  of g (dimension of the compact part of the Cartan decomposition for g).
##
DeclareGlobalFunction( "CompactDimension" );

#############################################################################
##
#F  Separate
##
##  The some technical function
##
DeclareGlobalFunction( "Separate" );

#############################################################################
##
#F  NicerSemisimpleSubalgebras(typeRank)
##
##  The function is based on LieAlgebraAndSubalgebras function from SLA Package.
##  It returns a list of types of semisimple subalgebras of a complex Lie algebra
##  (for a given type and rank). Remark: some types may correspond to several linearly
##  equivalence subalgebras.
##
DeclareGlobalFunction( "NicerSemisimpleSubalgebras" );

#############################################################################
##
#F  NumberOfSubalgebraClasses(typeRankG, typeRankH)
##
##  This function returns the number of linear equivalence classes for a given
##  type of subalgebra (which is the value of the second argument) in complex Lie
##  algebra (which is the value of the first argument).
##
DeclareGlobalFunction( "NumberOfSubalgebraClasses" );

#############################################################################
##
#F  PotentialSubalgebras("type",rank,id)
##
##  The input is a triple corresponding to a simple real Lie algebra g of rank 
## up to 8. This function returns a list of semisimple subalgebras in g as described
##  in Algorithm 1 (Bocheński, Jastrzębski, and Tralle 2019).
##
DeclareGlobalFunction( "PotentialSubalgebras" );

#############################################################################
##
#F  PotentialSubalgebraPairs("type",rank,id)
##
##  The input is a triple corresponding to a simple real Lie algebra g of rank
##  up to 8. This function returns subalgebra pairs in g as described in Algorithm 2
##  (Bocheński, Jastrzębski, and Tralle 2019).
##
DeclareGlobalFunction( "PotentialSubalgebraPairs" );

#############################################################################
##
#F  GetSymbolSimple(type, rank, id)
##
##  This function works up to rank 8 and returns string corresponding to the
##  symbol of real simple Lie algebra. Remark: we use notation from (Helgason 2001).
##
DeclareGlobalFunction( "GetSymbolSimple" );

#############################################################################
##
#F  GetSymbolSemisimple(tuple)
##
##  The input is a tuple corresponding to the semisimple real Lie algebra g
##  (represented in our notation). This function returns the symbol of g.
##
DeclareGlobalFunction( "GetSymbolSemisimple" );

#############################################################################
##
#F  GetSymbolPairs(output of PotentialSubalgebraPairs function)
##
##  This function returns symbols of the triples.
##
DeclareGlobalFunction( "GetSymbolPairs" );

#############################################################################
##
#F  RealFormByTuple(tuple)
##
##  The input is a tuple corresponding to a semisimple real Lie algebra g
##  (represented in our notation). This function returns a Lie algebra object
##  corresponding to g.
##
DeclareGlobalFunction( "RealFormByTuple" );

#############################################################################
##
#F  CheckTuple(tuple)
##
##  The inupt is a tuple. This function returns “true” when the tuple corresponds
##  to a real semisimple Lie algebra (it works only up to 8 simple components).
##
DeclareGlobalFunction( "CheckTuple" );


#E  ckforms.gd  . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here