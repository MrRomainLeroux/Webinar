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

assign("KM",1.2)
assign("Q",10)
assign("V",4)

MyModelEquations = ModelODEquations( list("Resp1" = expression( C1 ),
                                          "Resp2" = expression( C2 ) ) ,

                                     list("Deriv_C1" = expression( (Q/V2)*C2 - (Q/V1)*C1 - VMAX*(C1/V1)/(KM + (C1/V1)) - (CL/V1)*C1 ) ,
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
MyStatisticalModel = defineParameter( MyStatisticalModel, pCL )
MyStatisticalModel = defineParameter( MyStatisticalModel, pVMAX )
MyStatisticalModel = defineParameter( MyStatisticalModel, pV1 )

### Create and add the responses to the statistical model
MyStatisticalModel = addResponse( MyStatisticalModel, Response( "Resp1", Proportional( sigma_slope = 0.15 ) ) )
MyStatisticalModel = addResponse( MyStatisticalModel, Response( "Resp2", Proportional( sigma_slope = 0.15 ) ) )

### Finaly assign the statistical model to the project
MyProject = defineStatisticalModel( MyProject, MyStatisticalModel )

### Create a design
MyDesign= Design("Design")


#DOSE = c(0.03, 0.1, 0.3, 1, 3, 10)
# 6 groups 1/grp/dose

### For each arm create and add the sampling times for each response
Bras_test_group1 = Arm( name="Bras_test_group1", arm_size = 6 )

Bras_test_group1 = addSampling( Bras_test_group1, SamplingTimes( outcome = "Resp1",
                                                                 sample_time = c(c(1, 4)/24, 1, 3, 7, 14, 21) ) )

Bras_test_group1 = addSampling( Bras_test_group1, SamplingTimes( outcome = "Resp2",
                                                                 sample_time = c(c(1, 4)/24, 1, 3, 7, 14, 21) ) )

Bras_test_group1 = addAdministration( Bras_test_group1, Administration( outcome = "Resp1", time_dose = c(0), amount_dose = c(1000*0.1) ) )

Bras_test_group1 = setInitialConditions( Bras_test_group1, list( "C1" = expression( dose_Resp1 ), "C2" = 0 ) )

### Add the arm to the design
MyDesign = addArm( MyDesign, Bras_test_group1 )

### Add the design to the project
MyProject = addDesign( MyProject, MyDesign )

evaluationPop = EvaluatePopulationFIM( MyProject )











if(F){
plotPredictedResponses =

  ggplot( predictedResponses, aes(x = time   , y = Resp1 )) +

  geom_point() +

  geom_line() +

  scale_y_continuous(trans = 'log10', limits = c( 0.01, 1e4) )


print( plotPredictedResponses )
}



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










