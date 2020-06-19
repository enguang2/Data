# Produce some plots from previously-calculated quantites.

pdf("data1.pdf",onefile=TRUE)
par(mfrow=c(4,1),mar=c(2,4,4,2))

plot(df.us_all_national$Week_start,df.us_all_national$CoVNL63_pos.pct,xlab="",ylab="NL63"); title("positive percentage")
plot(df.us_all_national$Week_start,df.us_all_national$CoV229E_pos.pct,xlab="",ylab="229E")
plot(df.us_all_national$Week_start,df.us_all_national$CoVOC43_pos.pct,xlab="",ylab="OC43")
plot(df.us_all_national$Week_start,df.us_all_national$CoVHKU1_pos.pct,xlab="",ylab="HKU1")

plot(df.us_all_national$Week_start,df.us_all_national$CoVNL63_ili_x_pos_pct,xlab="",ylab="NL63"); title("ILI x positive percentage")
plot(df.us_all_national$Week_start,df.us_all_national$CoV229E_ili_x_pos_pct,xlab="",ylab="229E")
plot(df.us_all_national$Week_start,df.us_all_national$CoVOC43_ili_x_pos_pct,xlab="",ylab="OC43")
plot(df.us_all_national$Week_start,df.us_all_national$CoVHKU1_ili_x_pos_pct,xlab="",ylab="HKU1")

dev.off()

pdf("Reff1.pdf",onefile=TRUE)
par(mfrow=c(4,1),mar=c(2,4,4,2))

plot(CoV_ili_x_pos_pct_daily$CoVNL63,pch=20,xlab="",ylab="NL63"); title("interpolated data")
plot(CoV_ili_x_pos_pct_daily$CoV229E,pch=20,xlab="",ylab="229E")
plot(CoV_ili_x_pos_pct_daily$CoVOC43,pch=20,xlab="",ylab="OC43")
plot(CoV_ili_x_pos_pct_daily$CoVHKU1,pch=20,xlab="",ylab="HKU1")

plot(pmax(0.005,CoV_ili_x_pos_pct_daily$CoVNL63),type="l",xlab="",ylab="NL63",log="y"); title("interpolated data, log scale (min 0.005)")
plot(pmax(0.005,CoV_ili_x_pos_pct_daily$CoV229E),type="l",xlab="",ylab="229E",log="y")
plot(pmax(0.005,CoV_ili_x_pos_pct_daily$CoVOC43),type="l",xlab="",ylab="OC43",log="y")
plot(pmax(0.005,CoV_ili_x_pos_pct_daily$CoVHKU1),type="l",xlab="",ylab="HKU1",log="y")

plot(Reff.CoV_ili_x_pos_pct_SARS$Week_start,Reff.CoV_ili_x_pos_pct_SARS$CoVNL63,xlab="",ylab="NL63"); title("Smoothed weekly Reff, using SARS interval")
plot(Reff.CoV_ili_x_pos_pct_SARS$Week_start,Reff.CoV_ili_x_pos_pct_SARS$CoV229E,xlab="",ylab="229E")
plot(Reff.CoV_ili_x_pos_pct_SARS$Week_start,Reff.CoV_ili_x_pos_pct_SARS$CoVOC43,xlab="",ylab="OC43")
plot(Reff.CoV_ili_x_pos_pct_SARS$Week_start,Reff.CoV_ili_x_pos_pct_SARS$CoVHKU1,xlab="",ylab="HKU1")

plot(Reff.CoV_ili_x_pos_pct_SARS$Week_start,Reff.CoV_ili_x_pos_pct_SARS$CoVHKU1,xlab="",ylab="HKU1",col="red",type="l",
     ylim=c(0,3))
lines(Reff.CoV_ili_x_pos_pct_SARS$Week_start,Reff.CoV_ili_x_pos_pct_SARS$CoVOC43,xlab="",ylab="OC43",col="blue",type="l")
abline(h=c(0,1,2,3))
title("Reff for OC43 (blue) and HKU1 (red), using SARS interval")

dev.off()

pdf("model1.pdf",onefile=TRUE)

par(mfrow=c(4,1),mar=c(2,4,4,2))

HKU1 <- df.for_model.beta$strain=="CoVHKU1"
OC43 <- ! HKU1

plot(df.for_model.beta$Week_start[HKU1],df.for_model.beta$depletion.same_strain[HKU1],ylab="HKU1 same depletion",xlab="",col="red",pch=20,ylim=c(0,300))
title("Depletion of suseptables from same and opposite strain")
plot(df.for_model.beta$Week_start[HKU1],df.for_model.beta$depletion.opp_strain[HKU1],ylab="HKU1 opposite depletion",xlab="",col="red",pch=20,ylim=c(0,300))
plot(df.for_model.beta$Week_start[OC43],df.for_model.beta$depletion.same_strain[OC43],ylab="OC43 same depletion",xlab="",col="blue",pch=20,ylim=c(0,300))
plot(df.for_model.beta$Week_start[OC43],df.for_model.beta$depletion.opp_strain[OC43],ylab="OC43 opposite depletion",xlab="",col="blue",pch=20,ylim=c(0,300))

par(mfrow=c(4,1),mar=c(2,4,4,2))

pred.beta <- exp(predict(model.beta))
actual_HKU1 <- left_join(df.for_model.beta[HKU1,],Reff.CoV_ili_x_pos_pct_SARS,"Week_start")$CoVHKU1_R_ili_x_pos_pct_SARS
actual_OC43 <- left_join(df.for_model.beta[OC43,],Reff.CoV_ili_x_pos_pct_SARS,"Week_start")$CoVOC43_R_ili_x_pos_pct_SARS

plot(df.for_model.beta$Week_start[HKU1],pred.beta[HKU1],pch=19,col="red",ylab="HKU1",xlab="",ylim=c(0,4))
abline(h=0:4,col="gray")
points(df.for_model.beta$Week_start[HKU1],actual_HKU1,pch=20)
title("Reff and model predictions")
plot(df.for_model.beta$Week_start[OC43],pred.beta[OC43],pch=19,col="blue",ylab="OC43",xlab="",ylim=c(0,4))
abline(h=0:4,col="gray")
points(df.for_model.beta$Week_start[OC43],actual_OC43,pch=20)

plot(df.for_model.beta$Week_start[HKU1],actual_HKU1-pred.beta[HKU1],pch=20,col="red",ylab="HKU1",xlab="")
abline(h=0,col="gray")
title("Residuals")
plot(df.for_model.beta$Week_start[OC43],actual_OC43-pred.beta[OC43],pch=20,col="blue",ylab="HKU1",xlab="")
abline(h=0,col="gray")

par(mfrow=c(4,1))

cf_HKU1 <- coef(model.beta)[c("depletion.same_strain","depletion.opp_strain",
           "season2","season3","season4","season5")]
cf_OC43 <- cf_HKU1 + coef(model.beta)[c("strainCoVOC43:depletion.same_strain","strainCoVOC43:depletion.opp_strain",
           "season2:strainCoVOC43","season3:strainCoVOC43","season4:strainCoVOC43","season5:strainCoVOC43")]

plot (df.for_model.beta$Week_start[HKU1], exp(df.for_model.beta$depletion.same_strain[HKU1]*cf_HKU1["depletion.same_strain"]), xlab="",ylab="HKU1 same strain", pch=20,col="red", ylim=c(0.5,1))
abline(h=c(0.5,0.75,1.0)); title("Effects of depletions")
abline(v=sapply(c("2014-10-01","2015-04-30","2015-10-01","2016-04-30","2016-10-01","2017-04-30","2017-10-01","2018-04-30","2018-10-01","2019-04-30"),as.Date))
plot (df.for_model.beta$Week_start[HKU1], exp(df.for_model.beta$depletion.opp_strain[HKU1]*cf_HKU1["depletion.opp_strain"]), xlab="",ylab="HKU1 opp strain", pch=20,col="blue", ylim=c(0.5,1))
abline(h=c(0.5,0.75,1.0))
abline(v=sapply(c("2014-10-01","2015-04-30","2015-10-01","2016-04-30","2016-10-01","2017-04-30","2017-10-01","2018-04-30","2018-10-01","2019-04-30"),as.Date))
plot (df.for_model.beta$Week_start[OC43], exp(df.for_model.beta$depletion.same_strain[OC43]*cf_OC43["depletion.same_strain"]), xlab="",ylab="OC43 same strain", pch=20,col="blue",ylim=c(0.5,1))
abline(h=c(0.5,0.75,1.0))
abline(v=sapply(c("2014-10-01","2015-04-30","2015-10-01","2016-04-30","2016-10-01","2017-04-30","2017-10-01","2018-04-30","2018-10-01","2019-04-30"),as.Date))
plot (df.for_model.beta$Week_start[OC43], exp(df.for_model.beta$depletion.opp_strain[OC43]*cf_OC43["depletion.opp_strain"]), xlab="",ylab="OC43 opp strain", pch=20,col="red", ylim=c(0.5,1))
abline(h=c(0.5,0.75,1.0))
abline(v=sapply(c("2014-10-01","2015-04-30","2015-10-01","2016-04-30","2016-10-01","2017-04-30","2017-10-01","2018-04-30","2018-10-01","2019-04-30"),as.Date))

par(mfrow=c(4,1),mar=c(2,4,4,2))

seas_HKU1 <- c(0,cf_HKU1[-c(1,2)])
seff_HKU1 <- exp(log(pred.beta)[HKU1]-coef(model.beta)["(Intercept)"]-seas_HKU1[df.for_model.beta$season[HKU1]]
               -df.for_model.beta$depletion.same_strain[HKU1]*cf_HKU1["depletion.same_strain"]
               -df.for_model.beta$depletion.opp_strain[HKU1]*cf_HKU1["depletion.opp_strain"])
plot (df.for_model.beta$Week_start[HKU1], seff_HKU1, xlab="", ylab="HKU1", pch=20, ylim=c(1.0,1.5))
title("Seasonality effect (should be same for all years and strains)")
abline(h=c(1.0,1.25,1.5))
abline(v=sapply(c("2014-10-01","2015-04-30","2015-10-01","2016-04-30","2016-10-01","2017-04-30","2017-10-01","2018-04-30","2018-10-01","2019-04-30"),as.Date))

seas_OC43 <- c(0,cf_OC43[-c(1,2)])
seff_OC43 <- exp(log(pred.beta)[OC43]-coef(model.beta)["(Intercept)"]-coef(model.beta)["strainCoVOC43"]-seas_OC43[df.for_model.beta$season[OC43]]
               -df.for_model.beta$depletion.same_strain[OC43] * cf_OC43["depletion.same_strain"]
               -df.for_model.beta$depletion.opp_strain[OC43]  * cf_OC43["depletion.opp_strain"])
plot (df.for_model.beta$Week_start[OC43], seff_OC43, xlab="", ylab="OC43", pch=20, ylim=c(1.0,1.5))
abline(h=c(1.0,1.25,1.5))
abline(v=sapply(c("2014-10-01","2015-04-30","2015-10-01","2016-04-30","2016-10-01","2017-04-30","2017-10-01","2018-04-30","2018-10-01","2019-04-30"),as.Date))

dev.off()

pdf("effects-HKU1-15-16.pdf",height=8,width=8.5)
par(mar=c(2,4,4,2),xaxt="n")

sel <- 2
w <- df.for_model.beta$season[HKU1]==sel
plot (df.for_model.beta$season_week[HKU1][w], seff_HKU1[w], xlab="", ylab="Multiplicative effects on R", type="l", ylim=c(0.5,2.0), col="goldenrod", lwd=4, xaxp=c(1,30,6))
title("HKU1, weeks 1-33 of 2015-2016 season")
abline(h=c(0.5,0.75,1.0,1.25,1.5,1.75,2.0))
abline(v=seq(1,30,len=7))
abline(v=c(1,30),lwd=4)

points(1,exp(--seas_HKU1[sel]),pch=19)

lines (df.for_model.beta$season_week[HKU1][w], exp(df.for_model.beta$depletion.same_strain[HKU1]*cf_HKU1["depletion.same_strain"])[w], col="red", lwd=4)
lines (df.for_model.beta$season_week[HKU1][w], exp(df.for_model.beta$depletion.opp_strain[HKU1]*cf_HKU1["depletion.opp_strain"])[w], col="blue", lwd=4)

dev.off()

# Print model coefficients.

co <- coef(model.beta)
print(t(t(co)))

# Seasonal effects for OC43

print(t(t(co[c("season2","season3","season4","season5")] + co["strainCoVOC43"] +
          co[c("season2:strainCoVOC43","season3:strainCoVOC43","season4:strainCoVOC43","season5:strainCoVOC43")])))

# Depletion effects for OC43

print(t(t(co[c("depletion.same_strain","depletion.opp_strain")] +
          co[c("strainCoVOC43:depletion.same_strain","strainCoVOC43:depletion.opp_strain")])))
