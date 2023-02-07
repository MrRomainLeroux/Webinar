options(warn=-1)

# -----------------------------------------------
# Create PFIM project
# -----------------------------------------------

MyProject=PFIMProject(name = "PopED_Example02")

MyStatisticalModel = StatisticalModel()

# !! need to adjust .relStep for gradient computing
MyStatisticalModel = setParametersOdeSolver( MyStatisticalModel, list( atol=1e-16, rtol=1e-6, .relStep = 1e-12 ) )

# -----------------------------------------------
# model equations
# -----------------------------------------------

#ke = (CL/V1)
#k12 = (Q/V1)
#k21= (Q/V2)
#CP = (C1/V1)

MyModelEquations = ModelODEquations( list("RespPK" = expression( C1 ),
                                          "RespPD" = expression( C2 ) ) ,

                                     list("Deriv_C1" = expression( (Q/V2)*C2 - (Q/V1)*C1 - VMAX*(C1/V1)/(KM + (C1/V1)) - (CL/V1)*C1 ) ,
                                          "Deriv_C2" = expression( (Q/V1)*C1 -  (Q/V2)*C2 ) ) )

# -----------------------------------------------
# Assign the equations to the model
# -----------------------------------------------

MyStatisticalModel = defineModelEquations( MyStatisticalModel, MyModelEquations )

# -----------------------------------------------
# Define the variables of the ode model
# -----------------------------------------------

vC1 = ModelVariable( "C1/V1" )
vC2 = ModelVariable( "C2" )

MyStatisticalModel = defineVariable( MyStatisticalModel, vC1 )
MyStatisticalModel = defineVariable( MyStatisticalModel, vC2 )

# -----------------------------------------------
# Set mu and omega for each parameter
# -----------------------------------------------

pCL = ModelParameter( "CL", mu = 0.5,
                      omega = sqrt( 0.2 ),
                      distribution = LogNormalDistribution() )

pVMAX = ModelParameter( "VMAX", mu = 20,
                        omega = sqrt( 0.2 ),
                        distribution = LogNormalDistribution() )

pV1 = ModelParameter( "V1", mu = 2.5,
                      omega = sqrt( 0.1 ),
                      distribution = LogNormalDistribution() )

pQ = ModelParameter( "Q", mu = 10,
                      omega = sqrt( 0.0 ),
                      fixedMu = FALSE,
                      distribution = LogNormalDistribution() )

pV2 = ModelParameter( "V2", mu = 4,
                     omega = sqrt( 0.0 ),
                     fixedMu = FALSE,
                     distribution = LogNormalDistribution() )

pKM = ModelParameter( "KM", mu = 1.2,
                      omega = sqrt( 0.0 ),
                      fixedMu = FALSE,
                      distribution = LogNormalDistribution() )

# ------------------------------------------------
# Assign the parameters to the statistical model
# ------------------------------------------------

MyStatisticalModel = defineParameter( MyStatisticalModel, pCL )
MyStatisticalModel = defineParameter( MyStatisticalModel, pVMAX )
MyStatisticalModel = defineParameter( MyStatisticalModel, pV1 )
MyStatisticalModel = defineParameter( MyStatisticalModel, pQ )
MyStatisticalModel = defineParameter( MyStatisticalModel, pV2 )
MyStatisticalModel = defineParameter( MyStatisticalModel, pKM )

# -----------------------------------------------------
# Create and add the responses to the statistical model
# -----------------------------------------------------

MyStatisticalModel = addResponse( MyStatisticalModel, Response( "RespPK", Proportional( sigma_slope = 0.15 ) ) )

# -----------------------------------------------------
# Finaly assign the statistical model to the project
# -----------------------------------------------------

MyProject = defineStatisticalModel( MyProject, MyStatisticalModel )

# -----------------------------------------------------
# Create a design
# -----------------------------------------------------

MyDesign= Design("Design")

# -----------------------------------------------------------------------------------
# For each arm create and add the sampling times & administration for each response
# -----------------------------------------------------------------------------------

dose_group1 = 0.03
dose_group2 = 0.1
dose_group3 = 0.3
dose_group4 = 1
dose_group5 = 3
dose_group6 = 10

dose_group1 = 1000*dose_group1
dose_group2 = 1000*dose_group2
dose_group3 = 1000*dose_group3
dose_group4 = 1000*dose_group4
dose_group5 = 1000*dose_group5
dose_group6 = 1000*dose_group6

# sampleTimeGroup = c(c(1, 4)/24, 1, 3, 7, 14, 21)

# for plot original design
# !! system is computationally singular: reciprocal condition number = 1.1078e-25 for original design
# sampleTimeGroup =  c( c( (1:4)/24, 1:21 ) )

# group 1
Bras_test_group1 = Arm( name="Bras_test_group1", arm_size = 6 )
Bras_test_group1 = addSampling( Bras_test_group1, SamplingTimes( outcome = "RespPK",
                                                                 sample_time = sampleTimeGroup ) )
Bras_test_group1 = addSampling( Bras_test_group1, SamplingTimes( outcome = "RespPD",
                                                                 sample_time = sampleTimeGroup ) )
Bras_test_group1 = addAdministration( Bras_test_group1, Administration( outcome = "RespPK", time_dose = c(0), amount_dose = dose_group1 ) )
Bras_test_group1 = setInitialConditions( Bras_test_group1, list( "C1" = expression( dose_RespPK ), "C2" = 0 ) )

# group 2
Bras_test_group2 = Arm( name="Bras_test_group2", arm_size = 6 )
Bras_test_group2 = addSampling( Bras_test_group2, SamplingTimes( outcome = "RespPK",
                                                                 sample_time = sampleTimeGroup ) )
Bras_test_group2 = addSampling( Bras_test_group2, SamplingTimes( outcome = "RespPD",
                                                                 sample_time = sampleTimeGroup ) )
Bras_test_group2 = addAdministration( Bras_test_group2, Administration( outcome = "RespPK", time_dose = c(0), amount_dose = dose_group2 ) )
Bras_test_group2 = setInitialConditions( Bras_test_group2, list( "C1" = expression( dose_RespPK ), "C2" = 0 ) )

# group 3
Bras_test_group3 = Arm( name="Bras_test_group3", arm_size = 6 )
Bras_test_group3 = addSampling( Bras_test_group3, SamplingTimes( outcome = "RespPK",
                                                                 sample_time = sampleTimeGroup ) )
Bras_test_group3 = addSampling( Bras_test_group3, SamplingTimes( outcome = "RespPD",
                                                                 sample_time = sampleTimeGroup ) )
Bras_test_group3 = addAdministration( Bras_test_group3, Administration( outcome = "RespPK", time_dose = c(0), amount_dose = dose_group3 ) )
Bras_test_group3 = setInitialConditions( Bras_test_group3, list( "C1" = expression( dose_RespPK ), "C2" = 0 ) )

# group 4
Bras_test_group4 = Arm( name="Bras_test_group4", arm_size = 6 )
Bras_test_group4 = addSampling( Bras_test_group4, SamplingTimes( outcome = "RespPK",
                                                                 sample_time = sampleTimeGroup ) )
Bras_test_group4 = addSampling( Bras_test_group4, SamplingTimes( outcome = "RespPD",
                                                                 sample_time = sampleTimeGroup ) )
Bras_test_group4 = addAdministration( Bras_test_group4, Administration( outcome = "RespPK", time_dose = c(0), amount_dose = dose_group4 ) )
Bras_test_group4 = setInitialConditions( Bras_test_group4, list( "C1" = expression( dose_RespPK ), "C2" = 0 ) )

# group 5
Bras_test_group5 = Arm( name="Bras_test_group5", arm_size = 6 )
Bras_test_group5 = addSampling( Bras_test_group5, SamplingTimes( outcome = "RespPK",
                                                                 sample_time = sampleTimeGroup ) )
Bras_test_group5 = addSampling( Bras_test_group5, SamplingTimes( outcome = "RespPD",
                                                                 sample_time = sampleTimeGroup ) )
Bras_test_group5 = addAdministration( Bras_test_group5, Administration( outcome = "RespPK", time_dose = c(0), amount_dose = dose_group5 ) )
Bras_test_group5 = setInitialConditions( Bras_test_group5, list( "C1" = expression( dose_RespPK ), "C2" = 0 ) )

# group 6
Bras_test_group6 = Arm( name="Bras_test_group6", arm_size = 6 )
Bras_test_group6 = addSampling( Bras_test_group6, SamplingTimes( outcome = "RespPK",
                                                                 sample_time = sampleTimeGroup ) )
Bras_test_group6 = addSampling( Bras_test_group6, SamplingTimes( outcome = "RespPD",
                                                                 sample_time = sampleTimeGroup ) )
Bras_test_group6 = addAdministration( Bras_test_group6, Administration( outcome = "RespPK", time_dose = c(0), amount_dose = dose_group6 ) )
Bras_test_group6 = setInitialConditions( Bras_test_group6, list( "C1" = expression( dose_RespPK ), "C2" = 0 ) )

# Add the arm to the design
MyDesign = addArm( MyDesign, Bras_test_group1 )
MyDesign = addArm( MyDesign, Bras_test_group2 )
MyDesign = addArm( MyDesign, Bras_test_group3 )
MyDesign = addArm( MyDesign, Bras_test_group4 )
MyDesign = addArm( MyDesign, Bras_test_group5 )
MyDesign = addArm( MyDesign, Bras_test_group6 )

# Add the design to the project
MyProject = addDesign( MyProject, MyDesign )

# ---------------------------------------------------
# evaluate the population FIM
# ---------------------------------------------------

evaluationPop = EvaluatePopulationFIM( MyProject )

# ---------------------------------------------------
# display the results
# ---------------------------------------------------

show( evaluationPop )

# plot SE and RSE
plotSE(evaluationPop)

plotRSE(evaluationPop)

# !! plotResponse: logAxis unavalaible in the plot options
# plot evaluated design

designs = getDesign( evaluationPop )
nameDesign = getNameDesign(MyDesign)
resultsEvaluation = designs[[nameDesign]]@concentration
resultsResponse = do.call(rbind,resultsEvaluation)
lengthSamplingTimes = length( designs[[nameDesign]]@concentration$Bras_test_group1$time )
groupNames = paste0("Group",rep(1: getArmSize(Bras_test_group1) , each=lengthSamplingTimes))
resultsResponse = cbind( resultsResponse, groupNames )
colnames(resultsResponse) = c("time","RespPK","Group")

ggplot( resultsResponse, aes( x = time, y = RespPK , colour = Group ) ) +
  geom_point() +
  geom_line() +
  scale_x_continuous( "Time fro dose (days)" ) +
  scale_y_continuous( "Model predictions", trans = 'log10', limits = c( 0.01, 1e4) )





