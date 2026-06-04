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
only works with a Poisson observation model.

Note on the negative binomial dispersion parameter: 
Care is needed for the dispersion parameter, theta.d, prior. It can be weakly identified, particularly without abundant counts, 
and when not all marked individual samples are identified to individual (theta.marked[1] < 1). Identifiability is improved with 
telemetry and/or a marking process in generalized SMR. The model files are set up with uninformative priors for theta.d.
Another option currently commented out is a moderately informative prior for moderate to strong overdispersion. Using a single session
version in Spatial-Mark-Resight-Conditional, this prior worked well in one simulation scenario with negative binomial data
and did not produce bias in another simulation scenario where I simulated Poisson (theta.d very large) data. 
So, it did not appear very influential, but it could be with very sparse data.