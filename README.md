# Emergent-Compound-Learning

This is the repository for the Emergent Compound Learning, a final project for CCM.

## Contributors

Feng Cheng (fc1367@nyu.edu)

Hanxiao Lu (hl4860@nyu.edu)

Shannon Yasuda (sy3494@nyu.edu)

Shanshan Li (sl9712@nyu.edu)

## Structural Overview

### File related to the human experiment
- `Alchemy Setup.md` explains the structure of the alchemy game
- `Human.html`, `Human.js`, and the `/img/` directory, contain the interactive version of the alchemy game that the human subjects see. You can also try this game through `Human.html` if you have the three related items downloaded
- `survey_link` contains a link to the Qualtrics survey where the subjects write down their predicted products of all 11 compounds

### File related to data and results
- The `/Subject_Data/` directory contains processed and deidentified subject data
- The `/Subject_Analysis/` directory contains all the plots from the analysis
- The `/Figures/` directory contains the data and the plots from the model simulation

### File related to modeling
- The `/src/` directory contains all the code for modeling and data analysis
  - `Models.py` contains all the model classes
  - `human_analysis.py` contains the code for analyzing human data
