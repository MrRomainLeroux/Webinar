
# ******************************************************************

# Example 01 from PopED

# Evaluation

# ******************************************************************

# Create PFIM project
MyProject=PFIMProject(name = "PFIM_model_01")

# Create the statistical model
MyStatisticalModel=StatisticalModel()

# initial equations
#MyPKModel = ModelEquations( list( "RespPK1" = expression( ( dose_RespPK1 * ka/(V * (ka - Cl/V))) * (exp(-Cl/V *t) - exp(-ka * t ) ) ) ) )

# constant parameters
Wt = 32
WtCl = 0.75
WtV = 1

# change of variables
#Cl = (Cl*(Wt/70)**(WtCl))
#V  = (V*(Wt/70)**(WtV))

MyPKModel = ModelEquations( list( "RespPK1" = expression( ( dose_RespPK1 * ka/((V*(Wt/70)**(WtV)) *
                                                                                 (ka - (Cl*(Wt/70)**(WtCl))/(V*(Wt/70)**(WtV))))) *
                                                            (exp(-(Cl*(Wt/70)**(WtCl))/(V*(Wt/70)**(WtV)) *t) - exp(-ka * t ) ) ) ) )

# Assign the equations to the statistical model
MyStatisticalModel = defineModelEquations( MyStatisticalModel, MyPKModel)

# Set mu and omega for each parameter
pka = ModelParameter( "ka", mu = 0.25,
                      omega = sqrt( 0.08 ),
                      distribution = LogNormalDistribution() )

pV = ModelParameter( "V", mu = 100,
                     omega = sqrt( 0.1 ),
                     distribution = LogNormalDistribution() )

pCl = ModelParameter( "Cl", mu = 10,
                      omega = sqrt( 0.2 ),
                      distribution = LogNormalDistribution() )

# Assign the parameters to the statistical model
MyStatisticalModel = defineParameter( MyStatisticalModel, pka )
MyStatisticalModel = defineParameter( MyStatisticalModel, pV )
MyStatisticalModel = defineParameter( MyStatisticalModel, pCl )

# Create and add the responses to the statistical model
MyStatisticalModel = addResponse( MyStatisticalModel, Response( "RespPK1", Combined1( sigma_inter = 1.0, sigma_slope = 0.05 ) ) )

# Finaly assign the statistical model to the project
MyProject = defineStatisticalModel( MyProject, MyStatisticalModel )

# Create a design
MyDesign= Design("MyDesign")

# For each arm create and add the sampling times for each response
brasTest = Arm( name="Bras test", arm_size = 12 )

brasTest = addSampling( brasTest, SamplingTimes( outcome = "RespPK1",
                                                 sample_time = c( 5, 24, 48 ,72, 168 ) ) )

# sampling times c( 5, 24, 48, 72 )
# brasTest = addAdministration( brasTest, Administration( outcome = "RespPK1",
#                                                          time_dose = c( 5, 24, 48, 72 ),
#                                                          amount_dose = c( 100, 1000, 1000, 1000  ) ) )

# sampling times c( 5, 23, 47, 71 )
# amelioration des RSE
brasTest = addAdministration( brasTest, Administration( outcome = "RespPK1",
                                                        time_dose = c( 5, 23, 47, 71 ),
                                                        amount_dose = c( 100, 1000, 1000, 1000  ) ) )

# evaluation of the population FIM
evaluationPop = EvaluatePopulationFIM( MyProject )
show( evaluationPop )

# RSE
RSE = plotRSE( evaluationPop )
print( RSE )

# model response
plotResponse = plotResponse( evaluationPop, plotOptions = list( ) )
print( plotResponse )


##########################################################################################
# END example
##########################################################################################
