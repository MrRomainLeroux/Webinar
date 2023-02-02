
## https://metrumrg.com/course/a-gentle-introduction-to-optimaldesign-for-pharmacometric-models/


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
Wt = 70
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
MyStatisticalModel = addResponse( MyStatisticalModel, Response( "RespPK1",
                                                                Combined1( sigma_inter = 1.0, sigma_slope = 0.05 ) ) )

# Finaly assign the statistical model to the project
MyProject = defineStatisticalModel( MyProject, MyStatisticalModel )

# Create a design
MyDesign= Design("MyDesign")

# For each arm create and add the sampling times for each response
brasTest = Arm( name="Bras test", arm_size = 10 )

#brasTest = addSampling( brasTest, SamplingTimes( outcome = "RespPK1",
#                                                 sample_time = c( 1, 24, 48 , 96 ) ) )
# μ_ka              0.25  0.04531979  18.12792


#brasTest = addSampling( brasTest, SamplingTimes( outcome = "RespPK1",
#                                                 sample_time = c( 4, 24, 48 ,96 ) ) )
#μ_ka              0.25  0.05271620  21.08648


#brasTest = addSampling( brasTest, SamplingTimes( outcome = "RespPK1",
#                                                 sample_time = c( 6, 24, 48 ,96 ) ) )
#μ_ka              0.25  0.05919378  23.67751


#brasTest = addSampling( brasTest, SamplingTimes( outcome = "RespPK1",
#                                                sample_time = c( 12, 24, 48 ,96 ) ) )
#μ_ka              0.25  0.08969702  35.87881




brasTest = addAdministration( brasTest, Administration( outcome = "RespPK1",
                                                          time_dose = c( 0 ),
                                                          amount_dose = c( 10000 ) ) )

MyDesign <- addArm( MyDesign, brasTest )
### Add the design to the project
MyProject <- addDesign( MyProject, MyDesign )

# Evaluation of the Design
evaluationPop <- EvaluatePopulationFIM( MyProject )

show( evaluationPop )



# RSE
#RSE = plotRSE( evaluationPop )
#print( RSE )

# model response
plotResponse = plotResponse( evaluationPop, plotOptions = list( ) )
print( plotResponse )


##########################################################################################
# END example
##########################################################################################
