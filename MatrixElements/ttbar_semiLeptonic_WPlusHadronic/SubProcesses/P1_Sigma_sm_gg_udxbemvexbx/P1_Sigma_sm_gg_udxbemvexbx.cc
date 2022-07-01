// 
// *  This file was automatically generated by MoMEMta-MaGMEE,
// *  A MadGraph Matrix Element Exporter plugin for MoMEMta.
// *
// *  It is subject to MoMEMta-MaGMEE's license and copyright:
// *
// *  Copyright (C) 2016  Universite catholique de Louvain (UCL), Belgium
// *
// *  This program is free software: you can redistribute it and/or modify
// *  it under the terms of the GNU General Public License as published by
// *  the Free Software Foundation, either version 3 of the License, or
// *  (at your option) any later version.
// *
// *  This program is distributed in the hope that it will be useful,
// *  but WITHOUT ANY WARRANTY; without even the implied warranty of
// *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// *  GNU General Public License for more details.
// *
// *  You should have received a copy of the GNU General Public License
// *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
// 

#include <string> 
#include <utility> 
#include <vector> 
#include <map> 

#include <P1_Sigma_sm_gg_udxbemvexbx.h> 
#include <HelAmps_sm.h> 

#include <momemta/ParameterSet.h> 
#include <momemta/SLHAReader.h> 

namespace ttbar_semiLeptonic_WPlusHadronic_sm 
{

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g g > t t~ WEIGHTED<=2 @1
// *   Decay: t > w+ b WEIGHTED<=2
// *     Decay: w+ > u d~ WEIGHTED<=2
// *   Decay: t~ > w- b~ WEIGHTED<=2
// *     Decay: w- > e- ve~ WEIGHTED<=2
// Process: g g > t t~ WEIGHTED<=2 @1
// *   Decay: t > w+ b WEIGHTED<=2
// *     Decay: w+ > c s~ WEIGHTED<=2
// *   Decay: t~ > w- b~ WEIGHTED<=2
// *     Decay: w- > e- ve~ WEIGHTED<=2
// Process: u u~ > t t~ WEIGHTED<=2 @1
// *   Decay: t > w+ b WEIGHTED<=2
// *     Decay: w+ > u d~ WEIGHTED<=2
// *   Decay: t~ > w- b~ WEIGHTED<=2
// *     Decay: w- > e- ve~ WEIGHTED<=2
// Process: c c~ > t t~ WEIGHTED<=2 @1
// *   Decay: t > w+ b WEIGHTED<=2
// *     Decay: w+ > u d~ WEIGHTED<=2
// *   Decay: t~ > w- b~ WEIGHTED<=2
// *     Decay: w- > e- ve~ WEIGHTED<=2
// Process: d d~ > t t~ WEIGHTED<=2 @1
// *   Decay: t > w+ b WEIGHTED<=2
// *     Decay: w+ > u d~ WEIGHTED<=2
// *   Decay: t~ > w- b~ WEIGHTED<=2
// *     Decay: w- > e- ve~ WEIGHTED<=2
// Process: s s~ > t t~ WEIGHTED<=2 @1
// *   Decay: t > w+ b WEIGHTED<=2
// *     Decay: w+ > u d~ WEIGHTED<=2
// *   Decay: t~ > w- b~ WEIGHTED<=2
// *     Decay: w- > e- ve~ WEIGHTED<=2
// Process: u u~ > t t~ WEIGHTED<=2 @1
// *   Decay: t > w+ b WEIGHTED<=2
// *     Decay: w+ > c s~ WEIGHTED<=2
// *   Decay: t~ > w- b~ WEIGHTED<=2
// *     Decay: w- > e- ve~ WEIGHTED<=2
// Process: c c~ > t t~ WEIGHTED<=2 @1
// *   Decay: t > w+ b WEIGHTED<=2
// *     Decay: w+ > c s~ WEIGHTED<=2
// *   Decay: t~ > w- b~ WEIGHTED<=2
// *     Decay: w- > e- ve~ WEIGHTED<=2
// Process: d d~ > t t~ WEIGHTED<=2 @1
// *   Decay: t > w+ b WEIGHTED<=2
// *     Decay: w+ > c s~ WEIGHTED<=2
// *   Decay: t~ > w- b~ WEIGHTED<=2
// *     Decay: w- > e- ve~ WEIGHTED<=2
// Process: s s~ > t t~ WEIGHTED<=2 @1
// *   Decay: t > w+ b WEIGHTED<=2
// *     Decay: w+ > c s~ WEIGHTED<=2
// *   Decay: t~ > w- b~ WEIGHTED<=2
// *     Decay: w- > e- ve~ WEIGHTED<=2

//--------------------------------------------------------------------------

// Initialize process.

P1_Sigma_sm_gg_udxbemvexbx::P1_Sigma_sm_gg_udxbemvexbx(const ParameterSet&
    configuration)
{

  std::string param_card = configuration.get < std::string > ("card"); 
  params.reset(new Parameters_sm(SLHA::Reader(param_card))); 

  // Set external particle masses for this matrix element
  mME.push_back(std::ref(params->ZERO)); 
  mME.push_back(std::ref(params->ZERO)); 
  mME.push_back(std::ref(params->ZERO)); 
  mME.push_back(std::ref(params->ZERO)); 
  mME.push_back(std::ref(params->mdl_MB)); 
  mME.push_back(std::ref(params->ZERO)); 
  mME.push_back(std::ref(params->ZERO)); 
  mME.push_back(std::ref(params->mdl_MB)); 

  mapFinalStates[{2, -1, 5, 11, -12, -5}] = 
  {
    {
      &P1_Sigma_sm_gg_udxbemvexbx::matrix_1_gg_ttx_t_wpb_wp_udx_tx_wmbx_wm_emvex, 
      false, 
      {
        std::make_pair(21, 21), std::make_pair(21, 21)
      }, 
      256, 
      256
    }
    , 
    {
      &P1_Sigma_sm_gg_udxbemvexbx::matrix_1_uux_ttx_t_wpb_wp_udx_tx_wmbx_wm_emvex, 
      true, 
      {
        std::make_pair(2, -2), std::make_pair(4, -4), std::make_pair(1, -1),
            std::make_pair(3, -3), std::make_pair(2, -2), std::make_pair(4,
            -4), std::make_pair(1, -1), std::make_pair(3, -3)
      }, 
      256, 
      36
    }
  }; 

}

void P1_Sigma_sm_gg_udxbemvexbx::resetHelicities() 
{
  for (auto& finalState: mapFinalStates)
  {
    for (auto& subProcess: finalState.second)
    {
      subProcess.resetHelicities(); 
    }
  }
}


//--------------------------------------------------------------------------
// Evaluate |M|^2, return a map of final states

std::map < std::pair < int, int > , double >
    P1_Sigma_sm_gg_udxbemvexbx::compute(const std::pair < std::vector<double> ,
    std::vector<double> > &initialMomenta, const std::vector < std::pair < int,
    std::vector<double> > > &finalState)
{

  // Set initial particle momenta
  momenta[0] = (double * ) (&initialMomenta.first[0]); 
  momenta[1] = (double * ) (&initialMomenta.second[0]); 

  // Suppose final particles are passed in the "correct" order
  std::vector<int> selectedFinalState(8 - 2); 
  for (size_t index = 0; index < (8 - 2); index++ )
  {
    selectedFinalState[index] = finalState[index].first; 
    momenta[index + 2] = (double * ) (&finalState[index].second[0]); 
  }

  // Set the event specific parameters
  params->updateParameters(); 
  params->updateCouplings(); 

  // Initialise result object
  std::map < std::pair < int, int > , double > result; 

  // Define permutation
  int perm[8]; 
  for(int i = 0; i < 8; i++ )
  {
    perm[i] = i; 
  }

  for(auto &me: mapFinalStates[selectedFinalState])
  {

    double me_sum = 0; 
    double me_mirror_sum = 0; 

    for(int ihel = 0; ihel < 256; ihel++ )
    {

      if(me.goodHel[ihel])
      {

        double sum = 0.; 
        calculate_wavefunctions(perm, helicities[ihel]); 
        double meTemp = me.callback( * this); 
        sum += meTemp; 
        me_sum += meTemp/me.denominator; 

        if(me.hasMirrorProcess)
        {
          perm[0] = 1; 
          perm[1] = 0; 
          // Calculate wavefunctions
          calculate_wavefunctions(perm, helicities[ihel]); 
          // Mirror back
          perm[0] = 0; 
          perm[1] = 1; 
          meTemp = me.callback( * this); 
          sum += meTemp; 
          me_mirror_sum += meTemp/me.denominator; 
        }

        if( !sum)
          me.goodHel[ihel] = false; 
      }
    }

    for (auto const &initialState: me.initialStates)
    {
      result[initialState] = me_sum; 
      if (me.hasMirrorProcess)
        result[std::make_pair(initialState.second, initialState.first)] =
            me_mirror_sum;
    }
  }


  return result; 
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void P1_Sigma_sm_gg_udxbemvexbx::calculate_wavefunctions(const int perm[],
    const int hel[])
{
  // Calculate wavefunctions for all processes
  static std::complex<double> w[18][18]; 

  // Calculate all wavefunctions
  vxxxxx(&momenta[perm[0]][0], mME[0], hel[0], -1, w[0]); 
  vxxxxx(&momenta[perm[1]][0], mME[1], hel[1], -1, w[1]); 
  oxxxxx(&momenta[perm[2]][0], mME[2], hel[2], +1, w[2]); 
  ixxxxx(&momenta[perm[3]][0], mME[3], hel[3], -1, w[3]); 
  FFV2_3(w[3], w[2], params->GC_100, params->mdl_MW, params->mdl_WW, w[4]); 
  oxxxxx(&momenta[perm[4]][0], mME[4], hel[4], +1, w[5]); 
  FFV2_1(w[5], w[4], params->GC_100, params->mdl_MT, params->mdl_WT, w[6]); 
  oxxxxx(&momenta[perm[5]][0], mME[5], hel[5], +1, w[7]); 
  ixxxxx(&momenta[perm[6]][0], mME[6], hel[6], -1, w[8]); 
  FFV2_3(w[8], w[7], params->GC_100, params->mdl_MW, params->mdl_WW, w[9]); 
  ixxxxx(&momenta[perm[7]][0], mME[7], hel[7], -1, w[10]); 
  FFV2_2(w[10], w[9], params->GC_100, params->mdl_MT, params->mdl_WT, w[11]); 
  VVV1P0_1(w[0], w[1], params->GC_10, params->ZERO, params->ZERO, w[12]); 
  FFV1_1(w[6], w[0], params->GC_11, params->mdl_MT, params->mdl_WT, w[13]); 
  FFV1_2(w[11], w[0], params->GC_11, params->mdl_MT, params->mdl_WT, w[14]); 
  ixxxxx(&momenta[perm[0]][0], mME[0], hel[0], +1, w[15]); 
  oxxxxx(&momenta[perm[1]][0], mME[1], hel[1], -1, w[16]); 
  FFV1P0_3(w[15], w[16], params->GC_11, params->ZERO, params->ZERO, w[17]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV1_0(w[11], w[6], w[12], params->GC_11, amp[0]); 
  FFV1_0(w[11], w[13], w[1], params->GC_11, amp[1]); 
  FFV1_0(w[14], w[6], w[1], params->GC_11, amp[2]); 
  FFV1_0(w[11], w[6], w[17], params->GC_11, amp[3]); 

}
double P1_Sigma_sm_gg_udxbemvexbx::matrix_1_gg_ttx_t_wpb_wp_udx_tx_wmbx_wm_emvex() 
{

  static std::complex<double> ztemp; 
  static std::complex<double> jamp[2]; 
  // The color matrix
  static const double denom[2] = {1, 1}; 
  static const double cf[2][2] = {{16, -2}, {-2, 16}}; 

  // Calculate color flows
  static const std::complex<double> cI(0., 1.); 
  jamp[0] = +cI * amp[0] - amp[1]; 
  jamp[1] = -cI * amp[0] - amp[2]; 

  // Sum and square the color flows to get the matrix element
  double matrix = 0; 
  for(int i = 0; i < 2; i++ )
  {
    ztemp = 0.; 
    for(int j = 0; j < 2; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j]; 
    matrix = matrix + real(ztemp * conj(jamp[i]))/denom[i]; 
  }

  return matrix; 
}

double P1_Sigma_sm_gg_udxbemvexbx::matrix_1_uux_ttx_t_wpb_wp_udx_tx_wmbx_wm_emvex() 
{

  static std::complex<double> ztemp; 
  static std::complex<double> jamp[2]; 
  // The color matrix
  static const double denom[2] = {1, 1}; 
  static const double cf[2][2] = {{27, 9}, {9, 27}}; 

  // Calculate color flows
  jamp[0] = +1./2. * (+1./3. * amp[3]); 
  jamp[1] = +1./2. * (-amp[3]); 

  // Sum and square the color flows to get the matrix element
  double matrix = 0; 
  for(int i = 0; i < 2; i++ )
  {
    ztemp = 0.; 
    for(int j = 0; j < 2; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j]; 
    matrix = matrix + real(ztemp * conj(jamp[i]))/denom[i]; 
  }

  return matrix; 
}



}

// Register matrix element with MoMEMta
#include <momemta/MatrixElementFactory.h> 
REGISTER_MATRIX_ELEMENT("ttbar_semiLeptonic_WPlusHadronic_sm_P1_Sigma_sm_gg_udxbemvexbx", ttbar_semiLeptonic_WPlusHadronic_sm::P1_Sigma_sm_gg_udxbemvexbx); 

