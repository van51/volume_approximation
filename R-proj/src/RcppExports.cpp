// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// InnerBall
Rcpp::NumericVector InnerBall(Rcpp::NumericMatrix A, bool Zono, bool Vpoly);
RcppExport SEXP _volesti_InnerBall(SEXP ASEXP, SEXP ZonoSEXP, SEXP VpolySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< bool >::type Zono(ZonoSEXP);
    Rcpp::traits::input_parameter< bool >::type Vpoly(VpolySEXP);
    rcpp_result_gen = Rcpp::wrap(InnerBall(A, Zono, Vpoly));
    return rcpp_result_gen;
END_RCPP
}
// poly_gen
Rcpp::NumericMatrix poly_gen(int kind_gen, bool Vpoly_gen, int dim_gen, int m_gen);
RcppExport SEXP _volesti_poly_gen(SEXP kind_genSEXP, SEXP Vpoly_genSEXP, SEXP dim_genSEXP, SEXP m_genSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type kind_gen(kind_genSEXP);
    Rcpp::traits::input_parameter< bool >::type Vpoly_gen(Vpoly_genSEXP);
    Rcpp::traits::input_parameter< int >::type dim_gen(dim_genSEXP);
    Rcpp::traits::input_parameter< int >::type m_gen(m_genSEXP);
    rcpp_result_gen = Rcpp::wrap(poly_gen(kind_gen, Vpoly_gen, dim_gen, m_gen));
    return rcpp_result_gen;
END_RCPP
}
// rounding
Rcpp::NumericMatrix rounding(Rcpp::NumericMatrix A, unsigned int walk_len, bool coord, bool ball_walk, double delta, bool Vpoly, bool Zono);
RcppExport SEXP _volesti_rounding(SEXP ASEXP, SEXP walk_lenSEXP, SEXP coordSEXP, SEXP ball_walkSEXP, SEXP deltaSEXP, SEXP VpolySEXP, SEXP ZonoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< unsigned int >::type walk_len(walk_lenSEXP);
    Rcpp::traits::input_parameter< bool >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< bool >::type ball_walk(ball_walkSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< bool >::type Vpoly(VpolySEXP);
    Rcpp::traits::input_parameter< bool >::type Zono(ZonoSEXP);
    rcpp_result_gen = Rcpp::wrap(rounding(A, walk_len, coord, ball_walk, delta, Vpoly, Zono));
    return rcpp_result_gen;
END_RCPP
}
// Rsample_points
Rcpp::NumericMatrix Rsample_points(Rcpp::NumericMatrix A, unsigned int walk_len, double e, Rcpp::NumericVector InnerVec, bool CG, bool ball_walk, double delta, bool coord, bool Vpoly, bool Zono, bool sam_simplex, bool sam_can_simplex, bool sam_arb_simplex, bool sam_ball, bool sam_sphere, unsigned int numpoints, int dim_gen, double variance);
RcppExport SEXP _volesti_Rsample_points(SEXP ASEXP, SEXP walk_lenSEXP, SEXP eSEXP, SEXP InnerVecSEXP, SEXP CGSEXP, SEXP ball_walkSEXP, SEXP deltaSEXP, SEXP coordSEXP, SEXP VpolySEXP, SEXP ZonoSEXP, SEXP sam_simplexSEXP, SEXP sam_can_simplexSEXP, SEXP sam_arb_simplexSEXP, SEXP sam_ballSEXP, SEXP sam_sphereSEXP, SEXP numpointsSEXP, SEXP dim_genSEXP, SEXP varianceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< unsigned int >::type walk_len(walk_lenSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type InnerVec(InnerVecSEXP);
    Rcpp::traits::input_parameter< bool >::type CG(CGSEXP);
    Rcpp::traits::input_parameter< bool >::type ball_walk(ball_walkSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< bool >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< bool >::type Vpoly(VpolySEXP);
    Rcpp::traits::input_parameter< bool >::type Zono(ZonoSEXP);
    Rcpp::traits::input_parameter< bool >::type sam_simplex(sam_simplexSEXP);
    Rcpp::traits::input_parameter< bool >::type sam_can_simplex(sam_can_simplexSEXP);
    Rcpp::traits::input_parameter< bool >::type sam_arb_simplex(sam_arb_simplexSEXP);
    Rcpp::traits::input_parameter< bool >::type sam_ball(sam_ballSEXP);
    Rcpp::traits::input_parameter< bool >::type sam_sphere(sam_sphereSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type numpoints(numpointsSEXP);
    Rcpp::traits::input_parameter< int >::type dim_gen(dim_genSEXP);
    Rcpp::traits::input_parameter< double >::type variance(varianceSEXP);
    rcpp_result_gen = Rcpp::wrap(Rsample_points(A, walk_len, e, InnerVec, CG, ball_walk, delta, coord, Vpoly, Zono, sam_simplex, sam_can_simplex, sam_arb_simplex, sam_ball, sam_sphere, numpoints, dim_gen, variance));
    return rcpp_result_gen;
END_RCPP
}
// SliceSimplex
double SliceSimplex(Rcpp::NumericVector hyplane);
RcppExport SEXP _volesti_SliceSimplex(SEXP hyplaneSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type hyplane(hyplaneSEXP);
    rcpp_result_gen = Rcpp::wrap(SliceSimplex(hyplane));
    return rcpp_result_gen;
END_RCPP
}
// Rvolume
double Rvolume(Rcpp::NumericMatrix A, unsigned int walk_len, double e, Rcpp::NumericVector InnerVec, bool CG, unsigned int win_len, unsigned int N, double C, double ratio, double frac, bool ball_walk, double delta, bool Vpoly, bool Zono, bool coord, bool rounding);
RcppExport SEXP _volesti_Rvolume(SEXP ASEXP, SEXP walk_lenSEXP, SEXP eSEXP, SEXP InnerVecSEXP, SEXP CGSEXP, SEXP win_lenSEXP, SEXP NSEXP, SEXP CSEXP, SEXP ratioSEXP, SEXP fracSEXP, SEXP ball_walkSEXP, SEXP deltaSEXP, SEXP VpolySEXP, SEXP ZonoSEXP, SEXP coordSEXP, SEXP roundingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< unsigned int >::type walk_len(walk_lenSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type InnerVec(InnerVecSEXP);
    Rcpp::traits::input_parameter< bool >::type CG(CGSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type win_len(win_lenSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type C(CSEXP);
    Rcpp::traits::input_parameter< double >::type ratio(ratioSEXP);
    Rcpp::traits::input_parameter< double >::type frac(fracSEXP);
    Rcpp::traits::input_parameter< bool >::type ball_walk(ball_walkSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< bool >::type Vpoly(VpolySEXP);
    Rcpp::traits::input_parameter< bool >::type Zono(ZonoSEXP);
    Rcpp::traits::input_parameter< bool >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< bool >::type rounding(roundingSEXP);
    rcpp_result_gen = Rcpp::wrap(Rvolume(A, walk_len, e, InnerVec, CG, win_len, N, C, ratio, frac, ball_walk, delta, Vpoly, Zono, coord, rounding));
    return rcpp_result_gen;
END_RCPP
}
// vol_R
Rcpp::NumericMatrix vol_R(Rcpp::NumericMatrix A, unsigned int walk_len, double e, Rcpp::NumericVector InnerVec, bool CG, unsigned int win_len, unsigned int N, double C, double ratio, double frac, bool ball_walk, double delta, bool Vpoly, bool Zono, bool exact_zono, bool gen_only, bool Vpoly_gen, unsigned int kind_gen, unsigned int dim_gen, unsigned int m_gen, bool round_only, bool rotate_only, bool ball_only, bool sample_only, bool sam_simplex, bool sam_can_simplex, bool sam_arb_simplex, bool sam_ball, bool sam_sphere, unsigned int numpoints, double variance, bool construct_copula, Rcpp::NumericVector hyplane1, Rcpp::NumericVector hyplane2, unsigned int num_slices, bool sliceSimplex, bool coord, bool rounding, bool verbose);
RcppExport SEXP _volesti_vol_R(SEXP ASEXP, SEXP walk_lenSEXP, SEXP eSEXP, SEXP InnerVecSEXP, SEXP CGSEXP, SEXP win_lenSEXP, SEXP NSEXP, SEXP CSEXP, SEXP ratioSEXP, SEXP fracSEXP, SEXP ball_walkSEXP, SEXP deltaSEXP, SEXP VpolySEXP, SEXP ZonoSEXP, SEXP exact_zonoSEXP, SEXP gen_onlySEXP, SEXP Vpoly_genSEXP, SEXP kind_genSEXP, SEXP dim_genSEXP, SEXP m_genSEXP, SEXP round_onlySEXP, SEXP rotate_onlySEXP, SEXP ball_onlySEXP, SEXP sample_onlySEXP, SEXP sam_simplexSEXP, SEXP sam_can_simplexSEXP, SEXP sam_arb_simplexSEXP, SEXP sam_ballSEXP, SEXP sam_sphereSEXP, SEXP numpointsSEXP, SEXP varianceSEXP, SEXP construct_copulaSEXP, SEXP hyplane1SEXP, SEXP hyplane2SEXP, SEXP num_slicesSEXP, SEXP sliceSimplexSEXP, SEXP coordSEXP, SEXP roundingSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< unsigned int >::type walk_len(walk_lenSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type InnerVec(InnerVecSEXP);
    Rcpp::traits::input_parameter< bool >::type CG(CGSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type win_len(win_lenSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type C(CSEXP);
    Rcpp::traits::input_parameter< double >::type ratio(ratioSEXP);
    Rcpp::traits::input_parameter< double >::type frac(fracSEXP);
    Rcpp::traits::input_parameter< bool >::type ball_walk(ball_walkSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< bool >::type Vpoly(VpolySEXP);
    Rcpp::traits::input_parameter< bool >::type Zono(ZonoSEXP);
    Rcpp::traits::input_parameter< bool >::type exact_zono(exact_zonoSEXP);
    Rcpp::traits::input_parameter< bool >::type gen_only(gen_onlySEXP);
    Rcpp::traits::input_parameter< bool >::type Vpoly_gen(Vpoly_genSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type kind_gen(kind_genSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type dim_gen(dim_genSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type m_gen(m_genSEXP);
    Rcpp::traits::input_parameter< bool >::type round_only(round_onlySEXP);
    Rcpp::traits::input_parameter< bool >::type rotate_only(rotate_onlySEXP);
    Rcpp::traits::input_parameter< bool >::type ball_only(ball_onlySEXP);
    Rcpp::traits::input_parameter< bool >::type sample_only(sample_onlySEXP);
    Rcpp::traits::input_parameter< bool >::type sam_simplex(sam_simplexSEXP);
    Rcpp::traits::input_parameter< bool >::type sam_can_simplex(sam_can_simplexSEXP);
    Rcpp::traits::input_parameter< bool >::type sam_arb_simplex(sam_arb_simplexSEXP);
    Rcpp::traits::input_parameter< bool >::type sam_ball(sam_ballSEXP);
    Rcpp::traits::input_parameter< bool >::type sam_sphere(sam_sphereSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type numpoints(numpointsSEXP);
    Rcpp::traits::input_parameter< double >::type variance(varianceSEXP);
    Rcpp::traits::input_parameter< bool >::type construct_copula(construct_copulaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type hyplane1(hyplane1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type hyplane2(hyplane2SEXP);
    Rcpp::traits::input_parameter< unsigned int >::type num_slices(num_slicesSEXP);
    Rcpp::traits::input_parameter< bool >::type sliceSimplex(sliceSimplexSEXP);
    Rcpp::traits::input_parameter< bool >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< bool >::type rounding(roundingSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(vol_R(A, walk_len, e, InnerVec, CG, win_len, N, C, ratio, frac, ball_walk, delta, Vpoly, Zono, exact_zono, gen_only, Vpoly_gen, kind_gen, dim_gen, m_gen, round_only, rotate_only, ball_only, sample_only, sam_simplex, sam_can_simplex, sam_arb_simplex, sam_ball, sam_sphere, numpoints, variance, construct_copula, hyplane1, hyplane2, num_slices, sliceSimplex, coord, rounding, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_volesti_InnerBall", (DL_FUNC) &_volesti_InnerBall, 3},
    {"_volesti_poly_gen", (DL_FUNC) &_volesti_poly_gen, 4},
    {"_volesti_rounding", (DL_FUNC) &_volesti_rounding, 7},
    {"_volesti_Rsample_points", (DL_FUNC) &_volesti_Rsample_points, 18},
    {"_volesti_SliceSimplex", (DL_FUNC) &_volesti_SliceSimplex, 1},
    {"_volesti_Rvolume", (DL_FUNC) &_volesti_Rvolume, 16},
    {"_volesti_vol_R", (DL_FUNC) &_volesti_vol_R, 39},
    {NULL, NULL, 0}
};

RcppExport void R_init_volesti(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
