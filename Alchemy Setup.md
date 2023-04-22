## Alchemy Setup

### References:

- https://en.wikipedia.org/wiki/Alchemical_symbol
- https://en.wikipedia.org/wiki/List_of_alchemical_substances

### Background Story

In this setting, the subject is testing with a group of basic alchemical substances (elements) by throwing them into the universal solvant (*Alkahest*). The subject can select each element at most once per trial, and can throw in at most 4 elements per trial. This give us $2^4 = 16$ possible compounds, although only $15$ of them are meaningful because the empty compound (empty set) is just the solvent. 

Now, at each trial, the subject throw the elements into the solvent and wait for a bit for the reaction to occur. Once the reaction is settled, the subject perform a process called Virgo (Distillation) to get the product from the mixture. The reward value of each end product will be given to the subject, and will be noted in a notebook for later access. The subject can then move on to the next trial.

The goal of the experimentation is to identify the minimal configurations of all emergent compounds with the least amount of trials. We may add arbitrary costs to the simple elements and ask subject to maximize the reward (at which point it turns into a bandit problem), or simply ask subjects to deduce the emergent compound formulae.

### Basic Elements:

Below are the bask elements the subject sees.

| Element Name | Element Symbol |
| :----------: | :------------: |
|   Air (A)    |       🜁        |
|  Earth (E)   |       🜃        |
|   Fire (F)   |       🜂        |
|  Water (W)   |       🜄        |

### Emergent Compounds

Below are the product produced by emergent compounds.

| Compound Composition  |   Compound Product    | Product Symbol |
| :-------------------: | :-------------------: | :------------: |
|      $(\not 0)$       |         Null          |    $\not0$     |
|         $(A)$         | Antinomy (*Stibnium*) |       ♁        |
|         $(E)$         |  Bismuth (*Wismuth*)  |       ♆        |
|         $(F)$         | Arsenic (*Arsenicum*) |       🜺        |
|         $(W)$         |  Sulfur (*Sulphur*)   |       🜍        |
|     $(A\land E)$      |   Lead (*Plumbum*)    |       ♄        |
| $(A \land E \land W)$ |  Silver (*Argentum*)  |       ☽        |
| $(E \land F \land W)$ |    Gold (*Aurum*)     |       ☉        |

Note that, together with the unique compound, the yield of an emergent compound should also include the singular compound (so $(A \land E)$ produces Antinomy and Bismuth in addition to Lead). This is to keep the narrative consistent. That is, the singular compounds are not "used up" in creating novel compounds. Otherwise we will have a hard time explaining why a non-emergent compound produce singular compounds in addition to the composite compounds.
