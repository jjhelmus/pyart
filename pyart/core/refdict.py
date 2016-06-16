"""
pyart.core.refdict
==================

A dictionary-like class which references entries in another dictionary.

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    RefDict

"""

import collections
import itertools


class RefDict(collections.MutableMapping):
    """
    A dictionary-like class which references entries in another dictionary.

    Entries in this dictionary

    The comparison methods, __cmp__, __ge__, __gt__, __le__, __lt__, __ne__,
    nor the view methods, viewitems, viewkeys, viewvalues, are implemented.
    Neither is the the fromkeys method.

    Parameters
    ----------
    refdic : dict

    Examples
    --------
    >>> d = LazyLoadDict({'key1': 'value1', 'key2': 'value2'})
    >>> d.keys()
    ['key2', 'key1']
    >>> lazy_func = lambda : 999
    >>> d.set_lazy('lazykey1', lazy_func)
    >>> d.keys()
    ['key2', 'key1', 'lazykey1']
    >>> d['lazykey1']
    999

    """

    def __init__(self, refdic):
        """ initalize. """
        self._refdic = refdic
        self._nonref = {}

    # abstract methods
    def __setitem__(self, key, value):
        """ Set a entry. """
        if key in self._nonref:
            self._nonref[key] = value
        else:
            self._refdic[key] = value

    def __getitem__(self, key):
        """ Get the value of a key. """
        if key in self._nonref:
            return self._nonref[key]
        else:
            return self._refdic[key]

    def __delitem__(self, key):
        """ Remove a key from the dictionary. """
        if key in self._nonref:
            del self._nonref[key]
        else:
            del self._refdic[key]

    def __iter__(self):
        """ Iterate over all referenced and non-referenced keys. """
        return itertools.chain(self._refdic.copy(), self._nonref.copy())

    def __len__(self):
        """ Return the number keys. """
        all_keys = set(self._refdic.keys()).union(self._nonref.keys())
        return len(all_keys)

    # additional class to mimic dict behavior
    def __str__(self):
        """ Return a string representation of the object. """
        if len(self._refdic) == 0 or len(self._nonref) == 0:
            seperator = ''
        else:
            seperator = ', '
        lazy_strs = ['%s: LazyLoad(%s)' % r for r in lazy_reprs]
        lazy_str = ", ".join(lazy_strs) + '}'
        return str(self._refdic)[:-1] + seperator + lazy_str

    def has_key(self, key):
        """ True if dictionary has key, else False. """
        return key in self

    def copy(self):
        """
        Return a copy of the dictionary.

        Lazy keys are not evaluated in the original or copied dictionary.
        """
        dic = self.__class__(self._refdic.copy())
        # load all lazy keys into the copy
        for key, value_callable in self._nonref.items():
            dic.set_lazy(key, value_callable)
        return dic

    # lazy dictionary specific methods
    def set_nonref(self, key, value):
        """ Set a entry that is unique to this dictionary. """
        self._nonref[key] = value
