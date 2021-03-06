-- Available parameters : https://github.com/MoMEMta/MoMEMta/blob/master/core/src/MoMEMta.cc#L141
-- If "relative accuracy" is reached stop launching points unless "min_eval" point have been launched,
-- if after "max_eval" points have been launched stop launching even if we do not have reahced accuracy.
-- "n_start" : Vegas works by iterations will launch n_start point by n_start point and refine the grid at each step
-- n_increase allow to launch more points at each iteration
-- NB : all the point that have been launched are used to compute the intergal
cuba = {
    seed = 1,
    relative_accuracy = 0.01,
    verbosity = 1,
    max_eval = 10000000,
    n_start = 200000, 
    n_increase = 0,
}

-- NB: to be defined in the .cc is matrix_element_prefix

--inputs as will be feed to the blocks and ME
local neg_lepton = declare_input("lepton1")
local pos_lepton = declare_input("lepton2")
local bjet1 = declare_input("bjet1")
local bjet2 = declare_input("bjet2")
local isr = declare_input("isr")


TFFile = "/home/ucl/cp3/fbury/scratch/TransferFunctions/TF_norm.root"

parameters = {
    energy = 13000.,
    top_mass = 173.,
    top_width = 1.491500,
    W_mass = 80.419002,
    W_width = 2.047600,
    Z_mass = 91.1876,
    Z_width = 2.49,
    lep1_me_index = 1,
    lep1TFFile = TFFile,
    lep1TFName = "",
    lep2TFFile = TFFile,
    lep2TFName = "",
    jet1TFFile = TFFile,
    jet1TFName = "",
    jet2TFFile = TFFile,
    jet2TFName = "",

    export_graph_as = "dy_to_ll_simple_computing_graph.dot"
}

parameters.lep1TFName = {lep1_TF}
parameters.lep2TFName = {lep2_TF}
parameters.jet1TFName = {jet1_TF}
parameters.jet2TFName = {jet2_TF}

baseDirME = "/home/users/f/b/fbury/HHMEM/MEMWeight/MatrixElements/dy_to_llbb/"
matrix_element_lib = baseDirME .. "build/libme_dy_to_llbb.so"
ME_card = baseDirME .. "Cards/param_card.dat"
load_modules(matrix_element_lib)

-- Permute b and bbar
-- ME va etre moyenné sur les permutations --> effectuer les permutations uniquement si il y a une différence (sinon == *2/2) e.g. pas pour H--> bbbar
add_reco_permutations(bjet1, bjet2) 

-- Fix gen energy of the two leptons according to their transfer function
BinnedTransferFunctionOnEnergy.tf_neg_lepton = {
    ps_point = add_dimension(),
    reco_particle = neg_lepton.reco_p4,
    file = parameter('lep1TFFile'),
    th2_name = parameter('lep1TFName'),
}
neg_lepton.set_gen_p4("tf_neg_lepton::output")

BinnedTransferFunctionOnEnergy.tf_pos_lepton = {
    ps_point = add_dimension(),
    reco_particle = pos_lepton.reco_p4,
    file = parameter('lep2TFFile'),
    th2_name = parameter('lep2TFName'),
    --min_E = 5.,
}
pos_lepton.set_gen_p4("tf_pos_lepton::output")

-- Compute the phase space density for observed particles not concerned by the change of variable : here lepton on which we evaluate the TF 
-- The TF jacobian just account for dp/dE to go from |p| to E, not the phase space volume, the other blocks give the whole product of jacobian and phase space volume to go from one param to another
StandardPhaseSpace.phaseSpaceOut = {
    particles = {'tf_neg_lepton::output', 'tf_pos_lepton::output'}
}
-- Obtain the energy of the b's from momentum conservation (first two input are RECO b's, other inputs are used to fix PT to be balanced, all at gen level) 
BlockA.blocka = {
    p1 = bjet1.reco_p4,
    p2 = bjet2.reco_p4,
    branches = {'tf_neg_lepton::output', 'tf_pos_lepton::output', isr.reco_p4},
}
Looper.looper = {
    solutions = "blocka::solutions",
    path = Path("tfEval_bjet1", "tfEval_bjet2", "boost", "dy", "integrand") -- everything that will depend on the different solutions
}
bjet1.set_gen_p4("looper::particles/1")
bjet2.set_gen_p4("looper::particles/2")
BinnedTransferFunctionOnEnergyEvaluator.tfEval_bjet1 = {
    reco_particle = bjet1.reco_p4,
    gen_particle = bjet1.gen_p4, 
    file = parameter('jet1TFFile'),
    th2_name = parameter('jet1TFName'),
}
BinnedTransferFunctionOnEnergyEvaluator.tfEval_bjet2 = {
    reco_particle = bjet2.reco_p4,
    gen_particle = bjet2.gen_p4, 
    file = parameter('jet2TFFile'),
    th2_name = parameter('jet2TFName'),
}
genParticles = {
    pos_lepton.gen_p4,
    neg_lepton.gen_p4,
    bjet1.gen_p4,
    bjet2.gen_p4,
    }


BuildInitialState.boost = {
    do_transverse_boost = true,
    particles = genParticles
}
--P4Printer.printparton2 = {
--    input = 'boost::partons/2'
--}

jacobians = {'phaseSpaceOut::phase_space','tf_neg_lepton::TF_times_jacobian', 'looper::jacobian', 'tf_pos_lepton::TF_times_jacobian', 'tfEval_bjet1::TF', 'tfEval_bjet2::TF'}

MatrixElement.dy = {
  pdf = 'CT10nlo',
  pdf_scale = parameter('Z_mass'),

  matrix_element = 'dy_to_llbb_sm_P1_Sigma_sm_gg_epembbx',
  matrix_element_parameters = {
      card = ME_card
  },

  override_parameters = {
      mdl_MT = parameter('top_mass'),
  },

  initialState = 'boost::partons',

  particles = {
    inputs = genParticles,
    ids = {
      {
        pdg_id =  -11,
        me_index = 1,
      },

      {
        pdg_id = 11,
        me_index = 2,
      },

      {
        pdg_id = 5,
        me_index = 3,
      },

      {
        pdg_id = -5,
        me_index = 4,
      }
    }
  },

  jacobians = jacobians
}

DoubleLooperSummer.integrand = {
    input = "dy::output"
}

integrand("integrand::sum")
