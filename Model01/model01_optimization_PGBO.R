
# ******************************************************************

# Example 01 from PopED

# Design optimization

# ******************************************************************

### Create PFIM project
MyProject=PFIMProject(name = "PFIM_model_01")

### Create the statistical model
MyStatisticalModel=StatisticalModel()

Wt = 32
WtCl = 0.75
WtV = 1

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
                                                 sample_time = c( 5, 72, 96 ,120+168 ) ) )

brasTest = addAdministration( brasTest, Administration( outcome = "RespPK1",
                                                        time_dose = c( 0, 24, 48, 72, 96 ),
                                                        amount_dose = c( 10000, 10000, 10000, 10000,10000 ) ) )


# set optimization constraint

samplingBoundsConstraint = SamplingConstraint( response = "RespPK1", continuousSamplingTimes = list( c( 0, 6 ),
                                                                                                     c( 71, 72 ),
                                                                                                     c( 95, 96 ),
                                                                                                     c( 120 + 96, 120 + 168 ) ) )

Constr1 = DesignConstraint()
Constr1 = addSamplingConstraint( Constr1, samplingBoundsConstraint )
brasTest = addSamplingConstraints( brasTest, Constr1 )
MyProject = setConstraint( MyProject,Constr1 )

MyDesign = addArm( MyDesign, brasTest )
MyProject = addDesign( MyProject, MyDesign )

## run optimization
pgboOptimizer = PGBOAlgorithm( N = 200,
                               muteEffect = 0.5,
                               maxIteration = 100,
                               seed = 42,
                               showProcess = TRUE )

optimizationPGBOPopulationFIM = OptimizeDesign( MyProject, pgboOptimizer, PopulationFim() )

show( optimizationPGBOPopulationFIM )


# bpop[1] bpop[2] bpop[3] D[1,1] D[2,2] D[3,3] SIGMA[1,1] SIGMA[2,2]
# 13.38487 80.99149 119.36373 46.82132 253.41994 288.82097 24.18497 40.80823
#
# CL    = bpop[1] * exp(b[1]),
# V     = bpop[2] * exp(b[2]),
# KA    = bpop[3] * exp(b[3]),
# WT_CL = bpop[4],
# WT_V  = bpop[5],

