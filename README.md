# Reef Twin: the simulation model

Reef Twin is a reef aquarium simulator: a live tank you can keep, dose, stock and get wrong without risking real animals. This repository publishes the model behind it so anyone can check the work.

**Read the model: [MODEL.md](MODEL.md)**

It covers temperature and evaporation, carbonate chemistry (alkalinity, pH, calcium, magnesium), the nitrogen and phosphorus cycles, light, coral and fish biology, coralline algae, the cleanup crew, pests and equipment failures. Every constant is tagged with where it came from:

- **S**: a peer-reviewed or primary source
- **H**: hobby consensus with a named reference
- **E**: an honest estimate
- **T**: a tuning value

The claims that could not be sourced are listed together in §2.13, and the 18 scenarios the model must pass are in §2.12.

## Found something wrong?

Please [open an issue](https://github.com/michaeld812/reeftwin-model/issues). Tell me what the model does, what you'd expect instead, and your source or experience. Hobby experience counts, especially for the values tagged E.

## About

Written and maintained by Michael DeNigris, a software developer and reef keeper for over 20 years. The simulator's code is not public; this document is the model it implements, and it is updated when the model changes. See the commit history for what changed and when.
