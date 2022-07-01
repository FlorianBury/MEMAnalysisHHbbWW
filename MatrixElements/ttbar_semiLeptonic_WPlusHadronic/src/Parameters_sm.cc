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

#include <Parameters_sm.h> 

#include <momemta/SLHAReader.h> 

namespace ttbar_semiLeptonic_WPlusHadronic_sm 
{

using namespace std; 

// Constructor
Parameters_sm::Parameters_sm(const SLHA::Reader& card) 
{
// Prepare a vector for indices (for reading indexed entries)
vector<int> indices(2, 0); 

m_card_parameters["mdl_WH"] = card.get_block_entry("decay", 25, 6.382339e-03); 
m_card_parameters["mdl_WW"] = card.get_block_entry("decay", 24, 2.047600e+00); 
m_card_parameters["mdl_WZ"] = card.get_block_entry("decay", 23, 2.441404e+00); 
m_card_parameters["mdl_WT"] = card.get_block_entry("decay", 6, 1.491500e+00); 
m_card_parameters["mdl_ymtau"] = card.get_block_entry("yukawa", 15,
    1.777000e+00);
m_card_parameters["mdl_ymt"] = card.get_block_entry("yukawa", 6, 1.730000e+02); 
m_card_parameters["mdl_ymb"] = card.get_block_entry("yukawa", 5, 4.700000e+00); 
m_card_parameters["aS"] = card.get_block_entry("sminputs", 3, 1.180000e-01); 
m_card_parameters["mdl_Gf"] = card.get_block_entry("sminputs", 2,
    1.166390e-05);
m_card_parameters["aEWM1"] = card.get_block_entry("sminputs", 1, 1.325070e+02); 
m_card_parameters["mdl_MH"] = card.get_block_entry("mass", 25, 1.250000e+02); 
m_card_parameters["mdl_MZ"] = card.get_block_entry("mass", 23, 9.118800e+01); 
m_card_parameters["mdl_MTA"] = card.get_block_entry("mass", 15, 1.777000e+00); 
m_card_parameters["mdl_MT"] = card.get_block_entry("mass", 6, 1.730000e+02); 
m_card_parameters["mdl_MB"] = card.get_block_entry("mass", 5, 4.700000e+00); 

this->cacheParameters(); 
this->cacheCouplings(); 
}

void Parameters_sm::cacheParameters()
{
// Define "zero"
zero = 0; 
ZERO = 0; 
mdl_WH = m_card_parameters["mdl_WH"]; 
mdl_WW = m_card_parameters["mdl_WW"]; 
mdl_WZ = m_card_parameters["mdl_WZ"]; 
mdl_WT = m_card_parameters["mdl_WT"]; 
mdl_ymtau = m_card_parameters["mdl_ymtau"]; 
mdl_ymt = m_card_parameters["mdl_ymt"]; 
mdl_ymb = m_card_parameters["mdl_ymb"]; 
aS = m_card_parameters["aS"]; 
mdl_Gf = m_card_parameters["mdl_Gf"]; 
aEWM1 = m_card_parameters["aEWM1"]; 
mdl_MH = m_card_parameters["mdl_MH"]; 
mdl_MZ = m_card_parameters["mdl_MZ"]; 
mdl_MTA = m_card_parameters["mdl_MTA"]; 
mdl_MT = m_card_parameters["mdl_MT"]; 
mdl_MB = m_card_parameters["mdl_MB"]; mdl_conjg__CKM3x3 = 1.; 
mdl_CKM3x3 = 1.; 
mdl_conjg__CKM1x1 = 1.; 
mdl_complexi = std::complex<double> (0., 1.); 
mdl_MZ__exp__2 = ((mdl_MZ) * (mdl_MZ)); 
mdl_MZ__exp__4 = ((mdl_MZ) * (mdl_MZ) * (mdl_MZ) * (mdl_MZ)); 
mdl_sqrt__2 = sqrt(2.); 
mdl_MH__exp__2 = ((mdl_MH) * (mdl_MH)); 
mdl_aEW = 1./aEWM1; 
mdl_MW = sqrt(mdl_MZ__exp__2/2. + sqrt(mdl_MZ__exp__4/4. - (mdl_aEW * M_PI *
    mdl_MZ__exp__2)/(mdl_Gf * mdl_sqrt__2)));
mdl_sqrt__aEW = sqrt(mdl_aEW); 
mdl_ee = 2. * mdl_sqrt__aEW * sqrt(M_PI); 
mdl_MW__exp__2 = ((mdl_MW) * (mdl_MW)); 
mdl_sw2 = 1. - mdl_MW__exp__2/mdl_MZ__exp__2; 
mdl_cw = sqrt(1. - mdl_sw2); 
mdl_sqrt__sw2 = sqrt(mdl_sw2); 
mdl_sw = mdl_sqrt__sw2; 
mdl_g1 = mdl_ee/mdl_cw; 
mdl_gw = mdl_ee/mdl_sw; 
mdl_vev = (2. * mdl_MW * mdl_sw)/mdl_ee; 
mdl_vev__exp__2 = ((mdl_vev) * (mdl_vev)); 
mdl_lam = mdl_MH__exp__2/(2. * mdl_vev__exp__2); 
mdl_yb = (mdl_ymb * mdl_sqrt__2)/mdl_vev; 
mdl_yt = (mdl_ymt * mdl_sqrt__2)/mdl_vev; 
mdl_ytau = (mdl_ymtau * mdl_sqrt__2)/mdl_vev; 
mdl_muH = sqrt(mdl_lam * mdl_vev__exp__2); 
mdl_I1x33 = mdl_yb * mdl_conjg__CKM3x3; 
mdl_I2x33 = mdl_yt * mdl_conjg__CKM3x3; 
mdl_I3x33 = mdl_CKM3x3 * mdl_yt; 
mdl_I4x33 = mdl_CKM3x3 * mdl_yb; 
mdl_ee__exp__2 = ((mdl_ee) * (mdl_ee)); 
mdl_sw__exp__2 = ((mdl_sw) * (mdl_sw)); 
mdl_cw__exp__2 = ((mdl_cw) * (mdl_cw)); 
}

void Parameters_sm::cacheCouplings()
{
GC_100 = (mdl_ee * mdl_complexi * mdl_conjg__CKM1x1)/(mdl_sw * mdl_sqrt__2); 
}

void Parameters_sm::updateParameters()
{
mdl_sqrt__aS = sqrt(aS); 
G = 2. * mdl_sqrt__aS * sqrt(M_PI); 
mdl_G__exp__2 = ((G) * (G)); 
}

void Parameters_sm::updateCouplings()
{
GC_11 = mdl_complexi * G; 
GC_10 = -G; 
}

}

