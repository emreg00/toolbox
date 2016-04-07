
def keep_only_with_overlapping_values(key_to_values, set_to_overlap, n_min_overlap):
    """
    Keeping it non-one-liner for potential future feature additions based on values
    """
    key_to_values_mod = {}
    for key, values in key_to_values.iteritems():
	values &= set_to_overlap
	if len(values) >= n_min_overlap:
	    key_to_values_mod[key] = values
    return key_to_values_mod


def keep_only_with_overlapping_keys(key_to_values, dict_to_overlap):
    return dict((key, values) for key, values in key_to_values.iteritems() if key in dict_to_overlap)


def keep_only_with_unique_values(key_to_values, inner_delim="|"):
    values_to_keys = {}
    for key, values in key_to_values.iteritems():
	values_key = inner_delim.join(sorted(list(values)))
	values_to_keys.setdefault(values_key, []).append(key)
    key_to_values_mod = {}
    key_to_keys_equivalent = {}
    for values, keys in values_to_keys.iteritems():
	if len(keys) > 1:
	    keys.sort()
	    key_to_keys_equivalent[keys[0]] = set(keys[1:])
	key_to_values_mod[keys[0]] = set(values.split(inner_delim))
    return key_to_values_mod, key_to_keys_equivalent 


