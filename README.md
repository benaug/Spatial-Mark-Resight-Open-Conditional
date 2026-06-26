# Spatial-Mark-Resight-Open-Conditional

Jolly-Seber spatial mark-resight samplers that condition on and update the latent individual IDs. Integrated telemetry survival 

This model is an SMR version of the Jolly-Seber model here (JS-SCR-Dcov):

https://github.com/benaug/Jolly-Seber-N-Prior-DA

It also has allows integrated telemetry survival data, see here for an SCR version:

https://github.com/benaug/Spatial-IPM-Telemetry

The SMR model comes from the conditional SMR models here:

https://github.com/benaug/Spatial-Mark-Resight-Conditional

See here for the marginal version for open populations:

https://github.com/benaug/Spatial-Mark-Resight-Open-Marginal

I would use the marginal version unless there is overdispersion in the resighting counts. The marginal approach mixes better but
only works with a Poisson observation model. The mixing for the conditional version can be what you might call "terrible" in many
scenarios, so it may or may not be workable for any given data set, particularly when there is more posterior uncertainty.

There are 4 model versions with stationary activity centers--Poisson and negative binomial observation models and versions with and without interspersed
marking and sighting within years. There are 2 versions with RSF-based mobile activity centers--Poisson and negative binomial
observation models without interspersed marking and sighting. Mobile activity center test scripts are set up with abundant marked 
individuals with telemetry data, which is likely required for acceptable mixing.

Notes on the negative binomial dispersion parameter: 
Care is needed for the dispersion parameter, theta.d, prior. It can be weakly identified, particularly without abundant counts, 
and when not all marked individual samples are identified to individual (theta.marked[1] < 1). Identifiability is improved with 
telemetry and/or a marking process in generalized SMR. The NB model files are set up with moderately informative priors for moderate 
to strong overdispersion. This worked well in one simulation scenario and did not produce bias in another simulation scenario 
where I simulated Poisson (phi very large) data. So, it did not appear very influential, but it could be with very sparse data.
Can also try the commented out uniform prior.