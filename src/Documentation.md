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
	As such, when you put compounds into a query, make sure that they are in the correct order (in our case, `"A" < "B" < "C" < "D" < "E"`). If you are dealing with subject data (e.g., subject input), you might want to sort the tuple before giving it to the model like this 
	
	```python
	>>> test_tuple = ("C", "B", "A")
	>>> tuple(sorted(test))
	('A', 'B', 'C')
	```

In this implementation, the model of the environment is handled by `Models.Compound_Model`. This class keeps tracks of the subset/superset relationships between all compounds. It also handles reward calculation and allows the user to dynamically change the reward values of emergent compounds. As such, this class is intended to be used as a ground-truth model (or environment) and a interactive surface for some algorithmic model (say, a graph-searching algorithm).



## *Class* Models.Compound_Model

#### Class Attributes
> Note that the actual attributes are protected and the user only have access to the copied version. As such, later changes to the class attributes will not be reflected on these exported instances
- *elements*: A `np.ndarray` of `str`; contains elements
- *compounds*: A `np.ndarray` of `tuple`, contains compounds
- *E_compounds*: A `set` of `tuple`; contains emergent compounds
- *P_compounds*: A `set` of `tuple`; contains predicted compounds
- *subset_mat*: A  `np.ndarray` of `bool`; contains a square matrix storing pairwise subset relations (row as a subset of column)
- *superset_mat*: A  `np.ndarray` of `bool`; contains a square matrix storing pairwise superset relations (row as a superset of column)
- *compound_rewards*: A `dict` in the form of `tuple : float`; contain reward values of compounds

---

#### Compound_Model

`Models.Compound_Model(elements, emergent_compound_rewards = None)`

Initialize an instance.


***Input:***

- *elements*: A `np.ndarray` of `str`: contains elements
- *emergent_compound_rewards*:  A `dict` in the form of `tuple : float`; contain reward values of compounds. Can be omitted.

***output***:

- None

___


#### Compound_Model.add_emergent_compounds

`Compound_Model.add_emergent_compounds(emergent_compound_rewards)`

Add emergent compound rewards. Note that if an emergent compound is already defined in the instance, its old reward value will be overwritten. The instance will also dynamically update all compounds whose value are affected by this change.

***input:***

- *emergent_compound_rewards*:  A `dict` in the form of `tuple : float`; contain reward values of compounds.

***output***:

- None

---

#### Compound_Model.del_emergent_compounds

`Compound_Model.del_emergent_compounds(emergent_compounds)`:

Delete emergent compound rewards. The instance will also dynamically update all compounds whose value are affected by this change.

***input:***

- *emergent_compound*:  A `iterable` of `tuple`; contain emergent compounds to be deleted.

***output***:

- None

---

#### Compound_Model.initialize_rewards

`Compound_Model.initialize_rewards(emergent_compound_rewards)`:

Reinitializes all the reward values of compounds. Note that all class attributes related to reward (i.e.,`Compound_Model.E_compounds`, `Compound_Model.P_compounds`, `Compound_Model.compound_rewards`) are all reinitialized.

***input:***

- *emergent_compound_rewards*:  A `dict` in the form of `tuple : float`; contain reward values of compounds.

***output***:

- None

---

#### Compound_Model.get_rewards

`Compound_Model.get_rewards(query)`

Get the reward values of compounds contained in the query.

***input:***

- *query*:  An `tuple`, `int`, or an `iterable` containing `tuple` or `int` (`dict` is not accepted). When the innermost query item is a `tuple`, it must match a stored compound. If the innermost query item is an `int`, it must match the index of a stored compound (see **Compound_Model.compound_to_ind**).

***output***:

- A `np.ndarray` of `float`. If *query* is an `iterable`, the length of the output equals the length of the *query*.

---

#### Compound_Model.is_subset

`Compound_Model.is_subset(query_A, query_B)`

Get a matrix detailing whether compounds in query_A are subsets of compounds in query_B (whether the compound in the row is a subset of the compound in the column).

***input:***

- *query*:  An `tuple`, `int`, or an `iterable` containing `tuple` or `int` (`dict` is not accepted). When the innermost query item is a `tuple`, it must match a stored compound. If the innermost query item is an `int`, it must match the index of a stored compound (see **Compound_Model.compound_to_ind**).

***output***:

- A `np.ndarray` of `bool` of dim `len(query_A)*len(query_B)`

---

#### Compound_Model.is_superset

`Compound_Model.is_superset(query_A, query_B)`

Get a matrix detailing whether compounds in query_A are supersets of compounds in query_B (whether the compound in the row is a superset of the compound in the column).

***input:***

- *query*:  An `tuple`, `int`, or an `iterable` containing `tuple` or `int` (`dict` is not accepted). When the innermost query item is a `tuple`, it must match a stored compound. If the innermost query item is an `int`, it must match the index of a stored compound (see **Compound_Model.compound_to_ind**).

***output***:

- A `np.ndarray` of `bool` of dim `len(query_A)*len(query_B)`

---

#### Compound_Model.find_subsets

`Compound_Model.find_subsets(query, target = None, return_type = "c")`

Find all the compounds that are subsets of each compound in the query. This function can be used to decompose the reward contribution of a predicted compound to the rewards of emergent compounds by setting the emergent compounds as the target.

***input***:

-  *query*:  An `tuple`, `int`, or an `iterable` containing `tuple` or `int` (`dict` is not accepted). When the innermost query item is a `tuple`, it must match a stored compound. If the innermost query item is an `int`, it must match the index of a stored compound (see **Compound_Model.compound_to_ind**).
- *target*:  An `tuple`, `int`, or an `iterable` containing `tuple` or `int` (`dict` is not accepted). When the innermost query item is a `tuple`, it must match a stored compound. If the innermost query item is an `int`, it must match the index of a stored compound (see **Compound_Model.compound_to_ind**). If not defined, the function will set all the stored compounds as its target.
- *return_type*: a `str` in `("i", "c")`. Determine the return format as indices or compound tuples.

***output***:

- If `return_type == "c"`, return a `np.ndarray` of `tuple` as the compounds. Otherwise, return a `np.ndarray` of `int` as the indices of the compounds.

---

#### Compound_Model.find_supersets 

`Compound_Model.find_supersets(query, target = None, return_type = "c")`

Find all the compounds that are supersets of each compound in the query. This function can be used to find the compounds whose reward values are dependent on an emergent compound as the query.

***input***:

-  *query*:  An `tuple`, `int`, or an `iterable` containing `tuple` or `int` (`dict` is not accepted). When the innermost query item is a `tuple`, it must match a stored compound. If the innermost query item is an `int`, it must match the index of a stored compound (see **Compound_Model.compound_to_ind**).
- *target*:  An `tuple`, `int`, or an `iterable` containing `tuple` or `int` (`dict` is not accepted). When the innermost query item is a `tuple`, it must match a stored compound. If the innermost query item is an `int`, it must match the index of a stored compound (see **Compound_Model.compound_to_ind**). If not defined, the function will set all the stored compounds as its target.
- *return_type*: a `str` in `("i", "c")`. Determine the return format as indices or compound tuples.

***output***:

- If `return_type == "c"`, return a `np.ndarray` of `tuple` as the compounds. Otherwise, return a `np.ndarray` of `int` as the indices of the compounds.

---

#### Compound_Model.compound_to_ind

`Compound_Model.compound_to_ind(query)`

Find the indices of the compounds given their tuples.

***input***:

-  *query*:  An `tuple`, `int`, or an `iterable` containing `tuple` or `int` (`dict` is not accepted). When the innermost query item is a `tuple`, it must match a stored compound. If the innermost query item is an `int`, it must match the index of a stored compound (see **Compound_Model.compound_to_ind**).

***output***:

- A `np.ndarray` of `int` containing the indices

---

#### Compound_Model.ind_to_compount

`ind_to_compound(self, query)`

Find the tuples of the compounds given their indices.

***input***:

-  *query*:  An `tuple`, `int`, or an `iterable` containing `tuple` or `int` (`dict` is not accepted). When the innermost query item is a `tuple`, it must match a stored compound. If the innermost query item is an `int`, it must match the index of a stored compound (see **Compound_Model.compound_to_ind**).

***output***:

- A `np.ndarray` of `tuple` containing the compound tuples

---

#### Compound_Model.get_incidence_vectors

`get_incidence_vectors(self, query)`

Find the incidence vectors of the compounds. Incidence vectors represent whether an element is included in a compound. This form of representation can be convienent since many set operations can be acheved through linear algebra operations (which are much faster on numpy). 

 ***input***:

-  *query*:  An `tuple`, `int`, or an `iterable` containing `tuple` or `int` (`dict` is not accepted). When the innermost query item is a `tuple`, it must match a stored compound. If the innermost query item is an `int`, it must match the index of a stored compound (see **Compound_Model.compound_to_ind**).

***output***:

- A `np.ndarray` of `int` with dim of   `len(query)*len(Compound_Model.elements)`. `1` at pos `n` means that the element at pos `n` is present in the compound.

## Helper Functions

#### nparray_convert

`nparray_convert(arr)`:

A  function that helps to create an `np.ndarray` of object without letting numpy array constructor to convert the enclosed objects

***input***:

- *arr*: An `iterable`

***output***:

- A `np.ndarray`

---

#### multi_union

`multi_union(arr)`

A function that creates a set of union of members in the array. Useful when you want to find the union of several lists of compounds.

***input***:

- *arr*: An `iterable` containing `iterable` objects.

***output***:

- A `set`

---

#### subset_mat

`subset_mat(arr)`

A function that builds a pairwise subset relation matrix. The objects in `arr` must be hashable.

***input***:

- *arr*: An `iterable` containing hashable objects

***output***:

- A `np.ndarray` of dim `len(arr)*len(arr)`. The valus at `(i,j)` determines whether the object at pos `i` is a subset of the object at pos `j`

---

#### supserset_mat

`superset_mat(arr)`

A function that builds a pairwise superset relation matrix. The objects in `arr` must be hashable.

***input***:

- *arr*: An `iterable` containing hashable objects

***output***:

- A `np.ndarray` of dim `len(arr)*len(arr)`. The valus at `(i,j)` determines whether the object at pos `i` is a superset of the object at pos `j`

---

#### incidence_vectors

`incidence_vectors(elements, arr)`

A function that builds incidence vectors.

***input***:

- *elements*: An `iterable` containing `str`
- *arr*: An `iterable` containing `iterable`. 

***output***:

- An `np.ndarray` of dim `len(iterable)*len(elem)`. Each row vector at pos `i` determines whether the iterable at pos `i` in `arr` contains the elements.

---

#### powerset

`powerset(arr)`

A function that builds the powerset from a given iterable.

***input***:

- *arr*: An `iterable`

***output***:

- An `iterable` of `tuple` containing the powerset of the input `arr` (each tuple is a set in the powerset)

---

#### non_symmetrical_matrix_iteration

`non_symmetrical_matrix_iteration(data_array, target_matrix, function)`

A function wrapper that applys a given function to all pairs of elements in `data_array`. The functional output is non-symmetrical in the sense that $f(a,b) \not = f(b,a)$. This is useful for creating things like the subset matrix.

***input***:

- *data_array*: An `iterable`
- *target_matrix*: An `np.ndarray` of dim `len(data_array)*len(data_array)`. The result of the evaluation will be stored here.
- *function*: An `collable` object that takes two inputs and returns one output.

***output***: None

---

#### symmetrical_matrix_iteration

`symmetrical_matrix_iteration(data_array, target_matrix, function, skip_diagonal = True)`

A function wrapper that applys a given function to all pairs of elements in `data_array`. The functional output is symmetrical in the sense that $f(a,b) = f(b,a)$ (so only half of the matrix will need to be calculated). This is useful for creating things like a pairwise distance matrix.

***input***:

- *data_array*: An `iterable`
- *target_matrix*: An `np.ndarray` of dim `len(data_array)*len(data_array)`. The result of the evaluation will be stored here.
- *function*: An `collable` object that takes two inputs and returns one output.
- *skip_diagonal*: A `bool` that determines whether the function will be applied to calculate the values on the diagonal (a.k.a., $f(a,a)$). This is useful when the function cannot take two same inputs or when you want to define some arbitrary values for them.

***output***: None
