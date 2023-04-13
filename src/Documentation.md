# Documentation

This file contains simple documentation for  `Models.py`, which contains implementations of models and a dozen of helper functions.



## Overview

In the project, we treat elements as strings (e.g., `"A", "B", "C"`) and we treat compounds as *tuples* (e.g., `("A",), ("B", "C")`). This is because tuples are hashable and can be converted to dict keys. However, there are several things you should keep in mind when you work with compounds:

1. When you are creating a single-element compound (say, compound of "A"). You should create it in the form `("A",)` anod not simply `("A")`. Without a comma, python wil interpret `("A")` as `"A"`.
2. When you create an numpy array from a list of compound, you need to first define an empty numpy array with `dtype = object` and then individually move these compounds into the array. Simply converting a list of compounds by `np.array(compounds, dtype = object)` will end up with the compound tuples getting convereted into numpy arrays (which will be problematic because it's not hashable). You can use functions like `Model.nparray_convert` to do this.
3. Because tuples are ordered, two tuples with the same elements but in different order will not be considered identitcal: 
	```python
	>>> ("A", "B") == ("B", "A")
	False
	```
	As such, when you put compounds into a query, make sure that they are in the correct order (in our case, `"A" < "B" < "C" < "D" < "E"`).

In this implementation, the model of the environment is handled by `Models.Compound_Model`. This class keeps tracks of the subset/superset relationships between all compounds. It also handles reward calculation and allows the user to dynamically change the reward values of emergent compounds. As such, this class is intended to be used as a ground-truth model (or environment) and a interactive surface for some algorithmic model (say, a graph-searching algorithm).



## Classes
### Models.Compound_Model

#### Class Attributes
> Note that the actual attributes are protected and the user only have access to the copied version. As such, later changes to the class attributes will not be reflected on these exported instances
- *elements*: A `np.ndarray` of `str`; contains elements
- *compounds*: A `np.ndarray` of `tuple`, contains compounds
- *E_compounds*: A `set` of `tuple`; contains emergent compounds
- *P_compounds*: A `set` of `tuple`; contains predicted compounds
- *subset_mat*: A  `np.ndarray` of `bool`; contains a square matrix storing pairwise subset relations (row as a subset of column)
- *superset_mat*: A  `np.ndarray` of `bool`; contains a square matrix storing pairwise superset relations (row as a superset of column)
- *compound_rewards*: A `dict` in the form of `tuple : float`; contain reward values of compounds



#### Class Methods

**Models.Compound_Model** `(elements, emergent_compound_rewards = None)`:
Initialize an instance.
***Input:***

- *elements*: A `np.ndarray` of `str`: contains elements
- *emergent_compound_rewards*:  A `dict` in the form of `tuple : float`; contain reward values of compounds. Can be omitted.

***Output***:

- None



**Compound_Model.add_emergent_compounds** `(emergent_compound_rewards)`:

Add emergent compound rewards. Note that if an emergent compound is already defined in the instance, its old reward value will be overwritten. The instance will also dynamically update all compounds whose value are affected by this change.

***input:***

- *emergent_compound_rewards*:  A `dict` in the form of `tuple : float`; contain reward values of compounds.

***Output***:

- None



**Compound_Model.del_emergent_compounds**` (emergent_compounds)`:

Delete emergent compound rewards. The instance will also dynamically update all compounds whose value are affected by this change.

***input:***

- *emergent_compound*:  A `iterable` of `tuple`; contain emergent compounds to be deleted.

***Output***:

- None



**Compound_Model.initialize_rewards**` (emergent_compound_rewards)`:

Reinitializes all the reward values of compounds. Note that all class attributes related to reward (i.e.,`Compound_Model.E_compounds`, `Compound_Model.P_compounds`, `Compound_Model.compound_rewards`) are all reinitialized.

***input:***

- *emergent_compound_rewards*:  A `dict` in the form of `tuple : float`; contain reward values of compounds.

***Output***:

- None



**Compound_Model.get_rewards** `(query)`:

Get the reward values of compounds contained in the query.

***input:***

- *query*:  An `tuple`, `int`, or an `iterable` containing `tuple` or `int` (`dict` is not accepted). When the innermost query item is a `tuple`, it must match a stored compound. If the innermost query item is an `int`, it must match the index of a stored compound (see **Compound_Model.compound_to_ind**).

***Output***:

- A `np.ndarray` of `float`. If *query* is an `iterable`, the length of the output equals the length of the *query*.



**Compound_Model.is_subset** `(query_A, query_B)`:

Get a matrix detailing whether compounds in query_A are subsets of compounds in query_B (whether the compound in the row is a subset of the compound in the column).

***input:***

- *query*:  An `tuple`, `int`, or an `iterable` containing `tuple` or `int` (`dict` is not accepted). When the innermost query item is a `tuple`, it must match a stored compound. If the innermost query item is an `int`, it must match the index of a stored compound (see **Compound_Model.compound_to_ind**).

***Output***:

- A `np.ndarray` of `bool` of dim `len(query_A)*len(query_B)`



**Compound_Model.is_superset** `(query_A, query_B)`:

Get a matrix detailing whether compounds in query_A are supersets of compounds in query_B (whether the compound in the row is a superset of the compound in the column).

***input:***

- *query*:  An `tuple`, `int`, or an `iterable` containing `tuple` or `int` (`dict` is not accepted). When the innermost query item is a `tuple`, it must match a stored compound. If the innermost query item is an `int`, it must match the index of a stored compound (see **Compound_Model.compound_to_ind**).

***Output***:

- A `np.ndarray` of `bool` of dim `len(query_A)*len(query_B)`



**Compound_Model.find_subsets** `(query, target = None, return_type = "c")`:

Find all the compounds that are subsets of each compound in the query. This function can be used to decompose the reward contribution of a predicted compound to the rewards of emergent compounds by setting the emergent compounds as the target.

***input***:

-  *query*:  An `tuple`, `int`, or an `iterable` containing `tuple` or `int` (`dict` is not accepted). When the innermost query item is a `tuple`, it must match a stored compound. If the innermost query item is an `int`, it must match the index of a stored compound (see **Compound_Model.compound_to_ind**).
- *target*:  An `tuple`, `int`, or an `iterable` containing `tuple` or `int` (`dict` is not accepted). When the innermost query item is a `tuple`, it must match a stored compound. If the innermost query item is an `int`, it must match the index of a stored compound (see **Compound_Model.compound_to_ind**). If not defined, the function will set all the stored compounds as its target.
- *return_type*: a `str` in `("i", "c")`. Determine the return format as indices or compound tuples.

***output***:

- If `return_type == "c"`, return a `np.ndarray` of `tuple` as the compounds. Otherwise, return a `np.ndarray` of `int` as the indices of the compounds.



**Compound_Model.find_supersets** `(query, target = None, return_type = "c")`:

Find all the compounds that are supersets of each compound in the query. This function can be used to find the compounds whose reward values are dependent on an emergent compound as the query.

***input***:

-  *query*:  An `tuple`, `int`, or an `iterable` containing `tuple` or `int` (`dict` is not accepted). When the innermost query item is a `tuple`, it must match a stored compound. If the innermost query item is an `int`, it must match the index of a stored compound (see **Compound_Model.compound_to_ind**).
- *target*:  An `tuple`, `int`, or an `iterable` containing `tuple` or `int` (`dict` is not accepted). When the innermost query item is a `tuple`, it must match a stored compound. If the innermost query item is an `int`, it must match the index of a stored compound (see **Compound_Model.compound_to_ind**). If not defined, the function will set all the stored compounds as its target.
- *return_type*: a `str` in `("i", "c")`. Determine the return format as indices or compound tuples.

***output***:

- If `return_type == "c"`, return a `np.ndarray` of `tuple` as the compounds. Otherwise, return a `np.ndarray` of `int` as the indices of the compounds.
