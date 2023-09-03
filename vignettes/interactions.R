## ---- warning=FALSE, message=FALSE--------------------------------------------
library(ctmm) #load the package
data("buffalo") #import the data
projection(buffalo) <- median(buffalo) # reproject the data

# Fit the movement models to the tracking data
FITS <- list()
for(i in 1:length(buffalo))
{
  GUESS <- ctmm.guess(buffalo[[i]],interactive=FALSE)
  FITS[[i]] <- ctmm.select(buffalo[[i]],GUESS)
}
names(FITS) <- names(buffalo)

# calculate AKDES on a consistent grid
AKDES <- akde(buffalo,FITS,weights=TRUE,dt.plot=FALSE)

## -----------------------------------------------------------------------------
OVER <- overlap(AKDES)

## -----------------------------------------------------------------------------
class(OVER)
names(OVER)

## -----------------------------------------------------------------------------
OVER$CI

## -----------------------------------------------------------------------------
OVER$CI[,,"est"]

## -----------------------------------------------------------------------------
# pairwise CIs 
OVER$CI["Pepper","Toni",]
OVER$CI["Queen","Toni",]

## -----------------------------------------------------------------------------
plot(buffalo[c("Pepper", "Queen")],
     UD=AKDES[c("Pepper", "Queen")],
     col = c("#e76f51", "#264653"),
     col.DF=c("#f4a261", "#2a9d8f"),
     col.grid = NA)

## -----------------------------------------------------------------------------
OVER$CI["Pepper","Queen",]

## -----------------------------------------------------------------------------
CDE <- encounter(AKDES[c("Pepper", "Queen")])

#Visualise the CDE
plot(buffalo[c("Pepper", "Queen")],
     col=c("#e76f51", "#264653"),
     UD=CDE,
     col.DF="#046C9A",
     col.grid = NA)

## -----------------------------------------------------------------------------
plot(buffalo[c("Cilla", "Mvubu")],
     UD=AKDES[c("Cilla", "Mvubu")],
     col = c("#e76f51", "#264653"),
     col.DF=c("#f4a261", "#2a9d8f"),
     col.grid = NA)

## -----------------------------------------------------------------------------
DISTS <- distances(buffalo[c("Cilla","Mvubu")],FITS[c("Cilla","Mvubu")])

## -----------------------------------------------------------------------------
head(DISTS)

## -----------------------------------------------------------------------------
plot(DISTS$est ~ DISTS$timestamp,
     type = "l",
     col = "#5e548e",
     ylab = "Separation distance (m)",
     xlab = "")

## -----------------------------------------------------------------------------
cilla_sim <- simulate(FITS$Cilla,t=buffalo$Cilla$t)
mvubu_sim <- simulate(FITS$Mvubu,t=buffalo$Mvubu$t)

sim_dists <- distances(list(cilla_sim, mvubu_sim),FITS[c("Cilla","Mvubu")])

plot(list(cilla_sim, mvubu_sim),
     col = c("#e76f51", "#264653"),
     main = "Simulated data")

plot(sim_dists$est ~ sim_dists$timestamp,
     type = "l",
     col = "#5e548e",
     main = "Simulated distances",
     ylab = "Distance (m)",
     xlab = "Time",
     ylim = c(0,max(sim_dists$est)))

## -----------------------------------------------------------------------------
PROXIMITY <- proximity(buffalo[c("Cilla","Mvubu")],
                       FITS[c("Cilla","Mvubu")],
                       GUESS = ctmm(error=FALSE)) #this is to speed up the calculation

PROXIMITY

## -----------------------------------------------------------------------------
DISTS$encounter <- ifelse(DISTS$est <= 100, 1, 0)

## -----------------------------------------------------------------------------
plot(DISTS$encounter ~ DISTS$timestamp, xlab = "", ylab = "Encounter", main = "Scatter plot")

## -----------------------------------------------------------------------------
n <- sum(DISTS$encounter)
t <- "day" %#% (DISTS$t[nrow(DISTS)] - DISTS$t[1])
cat("There were an estimated ", n, " encounters between Cilla and Mvubu, and their encounter rate was ", round(n/t,2), " per day.")

## -----------------------------------------------------------------------------
enc_rad <- 1:1000
N <- vector("numeric", 1000)
for(i in 1:length(enc_rad)){
  N[i] <- sum(ifelse(DISTS$est <= enc_rad[i], 1, 0))
}

#visualise the results
plot(N ~ enc_rad,
     ylab = "Encounters",
     xlab = "Encounter radius (m)",
     type = "l",
     col = "#5e548e")

## -----------------------------------------------------------------------------
RATES <- rates(AKDES)

## -----------------------------------------------------------------------------
RATES$CI[,,"est"]

