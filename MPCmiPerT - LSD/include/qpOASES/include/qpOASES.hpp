
#if defined(__SINGLE_OBJECT__) || defined(__C_WRAPPER__)

#include <MessageHandling.cpp>
#include <Utils.cpp>
#include <Indexlist.cpp>
#include <SubjectTo.cpp>
#include <Bounds.cpp>
#include <Constraints.cpp>

#if !defined(__MATLAB__) || defined(WIN32)
#include <BLASReplacement.cpp>
#include <LAPACKReplacement.cpp>
#endif

#include <Matrices.cpp>
#include <Options.cpp>
#include <QProblemB.cpp>
#include <Flipper.cpp>
#include <QProblem.cpp>
#include <SQProblem.cpp>

#if defined(SOLVER_MA27) || defined(SOLVER_MA57)
#include <SparseSolver.cpp>
#include <SQProblemSchur.cpp>
#endif

#if !defined(__C_WRAPPER__) && !defined(__MATLAB__)
#include <OQPinterface.cpp>
#include <SolutionAnalysis.cpp>
#endif

#else /* default compilation mode */

#include <qpOASES/QProblemB.hpp>
#include <qpOASES/QProblem.hpp>
#include <qpOASES/SQProblem.hpp>
#include <qpOASES/SQProblemSchur.hpp>
#include <qpOASES/extras/OQPinterface.hpp>
#include <qpOASES/extras/SolutionAnalysis.hpp>

#endif
