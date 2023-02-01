
# ******************************************************************

# Example 01 from PopED

# Design optimization: PSO algorithm

# ******************************************************************

### Create PFIM project
MyProject=PFIMProject(name = "PFIM_model_01")

### Create the statistical model
MyStatisticalModel=StatisticalModel()

Wt = 32
WtCl = 0.75
WtV = 1

#Cl = (Cl*(Wt/70)**(WtCl))
#V  = (V*(Wt/70)**(WtV))

MyPKModel = ModelEquations( list( "RespPK1" = expression( ( dose_RespPK1 * ka/((V*(Wt/70)**(WtV)) *
                                                                                 (ka - (Cl*(Wt/70)**(WtCl))/(V*(Wt/70)**(WtV))))) *
                                                            (exp(-(Cl*(Wt/70)**(WtCl))/(V*(Wt/70)**(WtV)) *t) - exp(-ka * t ) ) ) ) )

### Assign the equations to the statistical model
MyStatisticalModel = defineModelEquations( MyStatisticalModel, MyPKModel)

### Set mu and omega for each parameter
pka = ModelParameter( "ka", mu = 0.25,
                      omega = sqrt( 0.08 ),
                      distribution = LogNormalDistribution() )

pV = ModelParameter( "V", mu = 100,
                     omega = sqrt( 0.1 ),
                     distribution = LogNormalDistribution() )

pCl = ModelParameter( "Cl", mu = 10,
                      omega = sqrt( 0.2 ),
                      distribution = LogNormalDistribution() )

### Assign the parameters to the statistical model
MyStatisticalModel = defineParameter( MyStatisticalModel, pka )
MyStatisticalModel = defineParameter( MyStatisticalModel, pV )
MyStatisticalModel = defineParameter( MyStatisticalModel, pCl )

### Create and add the responses to the statistical model
MyStatisticalModel = addResponse( MyStatisticalModel, Response( "RespPK1", Combined1( sigma_inter = 1.0, sigma_slope = 0.05 ) ) )

### Finaly assign the statistical model to the project
MyProject = defineStatisticalModel( MyProject, MyStatisticalModel )

### Create a design
MyDesign= Design("MyDesign")

### For each arm create and add the sampling times for each response
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

# set optimization constraint
samplingBoundsConstraint = SamplingConstraint( response = "RespPK1", continuousSamplingTimes = list( c( 0, 6 ),
                                                                                                     c( 23, 24 ),
                                                                                                     c( 46, 48 ),
                                                                                                     c( 69, 72 ),
                                                                                                     c( 96, 168 ) ) )
# add min_delay constraint
# no min_delay constraint by default !
samplingMinimalDelayConstraintRespPK1 = SamplingConstraint( response = "RespPK1", min_delay = 1 )

Constr1 = DesignConstraint()
Constr1 = addSamplingConstraint( Constr1, samplingBoundsConstraint )
Constr1 = addSamplingConstraint( Constr1, samplingMinimalDelayConstraintRespPK1 )

brasTest<- addSamplingConstraints( brasTest, Constr1 )

MyProject = setConstraint( MyProject,Constr1 )

MyDesign = addArm( MyDesign, brasTest )
MyProject = addDesign( MyProject, MyDesign )

## run optimization

# PSO algorithm
psoOptimizer = PSOAlgorithm( maxIteration = 500,
                             populationSize = 100,
                             personalLearningCoefficient = 2.05,
                             globalLearningCoefficient = 2.05,
                             showProcess = TRUE )

optimizationPSOPopulationFIM = OptimizeDesign( MyProject, psoOptimizer, PopulationFim() )

show( optimizationPSOPopulationFIM )

##########################################################################################
# END example
##########################################################################################
