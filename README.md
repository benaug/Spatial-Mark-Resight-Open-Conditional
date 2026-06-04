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