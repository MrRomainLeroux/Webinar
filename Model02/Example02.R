options(warn=-1)

### Create PFIM project
MyProject=PFIMProject(name = "PopED_Example02")

MyStatisticalModel = StatisticalModel()

#MyStatisticalModel = setParametersOdeSolver( MyStatisticalModel, list( .relStep = 1e-9 ) )

### model equations

#ke = (CL/V1)
#k12 = (Q/V1)
#k21= (Q/V2)
#CP = (C1/V1)

KM = 1.2
Q = 10
V2 = 4

MyModelEquations = ModelODEquations( list("Resp1" = expression( C1 ),
                                           "Resp2" = expression( C2 ) ) ,

                                      list("Deriv_C1" = expression( (Q/V2)*C2 - (Q/V1)*C1 - VMAX*CP/(KM + (C1/V1)) - (CL/V1)*C1 ) ,
                                           "Deriv_C2" = expression( (Q/V1)*C1 -  (Q/V2)*C2 ) ) )

### Assign the equations to the model
MyStatisticalModel = defineModelEquations( MyStatisticalModel, MyModelEquations )

### Define the variables of the ode model
vC1 = ModelVariable( "C1" )
vC2 = ModelVariable( "C2" )

MyStatisticalModel = defineVariable( MyStatisticalModel, vC1 )
MyStatisticalModel = defineVariable( MyStatisticalModel, vC2 )

#### Set mu and omega for each parameter

pCL = ModelParameter( "CL", mu = 0.5,
                      omega = sqrt( 0.2 ),
                      distribution = LogNormalDistribution() )

pVMAX = ModelParameter( "VMAX", mu = 20,
                        omega = sqrt( 0.2 ),
                        distribution = LogNormalDistribution() )

pV1 = ModelParameter( "V1", mu = 2.5,
                      omega = sqrt( 0.1 ),
                      distribution = LogNormalDistribution() )

### Assign the parameters to the statistical model
MyStatisticalModel = defineParameter( MyStatisticalModel, pV1 )
MyStatisticalModel = defineParameter( MyStatisticalModel, pRin )
MyStatisticalModel = defineParameter( MyStatisticalModel, pVMAX )



### Create and add the responses to the statistical model
MyStatisticalModel = addResponse( MyStatisticalModel, Response( "Resp1", Combined1( sigma_inter = 0.6, sigma_slope = 0.07 ) ) )
MyStatisticalModel = addResponse( MyStatisticalModel, Response( "Resp2", Proportionnal( sigma_slope = 4 ) ) )

### Finaly assign the statistical model to the project
MyProject = defineStatisticalModel( MyProject, MyStatisticalModel )

### Create a design
MyDesign= Design("MyDesign")

### For each arm create and add the sampling times for each response
brasTest = Arm( name="Bras test", arm_size = 32 )

brasTest = addSampling( brasTest, SamplingTimes( outcome = "Resp1",
                                                  sample_time =  c( 0.5, 1, 2, 6, 9, 12, 24, 36, 48, 72, 96, 120 ) ) )

brasTest = addSampling( brasTest, SamplingTimes( outcome = "Resp2",
                                                  sample_time = c( 0,24,36,48,72,96,120 ) ) )

brasTest = addAdministration( brasTest, Administration( outcome = "Resp1", time_dose = c(0), amount_dose = c(100) ) )

brasTest = setInitialConditions( brasTest, list( "C1" = expression( dose_RespPK/V ),
                                                  "C2" = expression( Rin/kout ) ) )

### Add the arm to the design
MyDesign = addArm( MyDesign, brasTest )

### Add the design to the project
MyProject = addDesign( MyProject, MyDesign )

evaluationPop = EvaluatePopulationFIM( MyProject )













if(F){
evaluationInd = EvaluateIndividualFIM( MyProject )
evaluationBay = EvaluateBayesianFIM( MyProject )

### Summary

show( evaluationPop )
show( evaluationInd )
show( evaluationBay )

outputPath = "C:/Users/ADMIN Romain LEROUX/Documents/GIT PFIM/PFIM/PFIM_CRAN/Tests_package/evaluation/ode"

plotOptions = list( unitTime=c("unit time"),
                    unitResponses= c("unit RespPK","unit RespPD" ) )

reportPFIMProject( evaluationPop,
                   outputPath = outputPath, plotOptions = plotOptions )

evaluationInd = setNamePFIMProject( evaluationInd, "PKPD_ode_bolus_individualFIM" )
reportPFIMProject( evaluationInd,
                   outputPath = outputPath, plotOptions = plotOptions )

evaluationBay = setNamePFIMProject( evaluationBay, "PKPD_ode_bolus_bayesianFIM" )
reportPFIMProject( evaluationBay,
                   outputPath = outputPath, plotOptions = plotOptions )
}










