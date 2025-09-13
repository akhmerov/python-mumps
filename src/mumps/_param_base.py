"""Private base classes for MUMPS parameter arrays and enum-like access.

This module provides:
- BoundToken: value objects representing a specific allowed choice; support .set(), int(), and docs
- BoundParam: descriptor-bound view that exposes tokens and current value
- ParamDescriptor and @param decorator: to define per-parameter metadata and tokens
- _OneBasedArray: minimal 1-based indexable array wrapper
- ParamArray: base class for ICNTL/CNTL-like arrays with validation via ParamDescriptor
- RawArray: simple 1-based storage for float/int arrays without validation
"""

from __future__ import annotations

from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple


class BoundToken:
    """A token bound to a specific parameter instance; supports int(), set(), and docs."""

    def __init__(
        self,
        owner_alias: str,
        param_name: str,
        code: int,
        doc: str,
        setter: Callable[[int], None],
    ):
        self._owner_alias = owner_alias  # e.g., "icntl"
        self._param_name = param_name  # e.g., "message_level"
        self._code = int(code)
        self._setter = setter
        self.__doc__ = doc

    def set(self) -> None:
        """Set the owning parameter to this token's code."""
        self._setter(self._code)

    def __int__(self) -> int:  # for assertions
        return self._code

    def __repr__(self) -> str:
        return (
            f"<{self._owner_alias}.{self._param_name}.{self._token_name}: {self._code}>"
        )

    # Internal: set by the factory so __repr__ shows the token's attribute name
    _token_name: str = ""


class BoundParam:
    """A bound parameter descriptor that exposes tokens as attributes and reflects current value.

    Provides:
    - int() -> current numeric code
    - .__doc__ -> concise param documentation
    - attribute access for tokens, each with .set() and docs
    """

    def __init__(
        self,
        *,
        owner_alias: str,
        param_name: str,
        idx1: int,
        get_code: Callable[[], float | int],
        set_code: Callable[[float | int], None],
        doc: str,
        tokens: Iterable[Tuple[str, int, str]],
        page: Optional[int] = None,
    ):
        self._owner_alias = owner_alias
        self._param_name = param_name
        self._idx1 = int(idx1)
        self._get_code = get_code
        self._set_code = set_code
        # Build docstring with optional page reference
        if page is not None:
            page_line = f"(See MUMPS Users' Guide, page {int(page)})"
            base = (doc or "").strip()
            self.__doc__ = (page_line + ("\n" + base if base else "")).strip()
        else:
            self.__doc__ = doc
        # Create stable token objects as attributes
        self._tokens: Dict[str, BoundToken] = {}
        for name, code, tdoc in tokens:
            tok = BoundToken(
                owner_alias, param_name, int(code), tdoc, setter=self._set_code
            )
            tok._token_name = name
            self._tokens[name] = tok
            setattr(self, name, tok)

    def __int__(self) -> int:
        return int(self._get_code())

    def __float__(self) -> float:
        val = self._get_code()
        try:
            return float(val)
        except Exception:
            return float(int(val))

    def _current_token_name(self) -> Optional[str]:
        code = int(self)
        for name, tok in self._tokens.items():
            if int(tok) == code:
                return name
        return None

    def __repr__(self) -> str:
        # Enumerated parameter: show token name when available, else numeric code
        if self._tokens:
            name = self._current_token_name()
            code = int(self)
            if name is None:
                return f"<{self._owner_alias}.{self._param_name}: {code}>"
            return f"<{self._owner_alias}.{self._param_name}.{name}: {code}>"
        # Non-enumerated: show raw value (float or int)
        val = self._get_code()
        return f"<{self._owner_alias}.{self._param_name}: {val}>"

    def __dir__(self) -> List[str]:
        # Expose token names for tab-completion
        base = set(super().__dir__())
        return sorted(base.union(self._tokens.keys()))

    def __getattr__(self, name: str) -> Any:
        # Provide dynamic access to tokens as attributes
        tok = self._tokens.get(name)
        if tok is not None:
            return tok
        raise AttributeError(name)


class ParamDescriptor:
    """Data descriptor representing a single parameter ICNTL(k), etc."""

    def __init__(
        self,
        *,
        index: int,
        doc: str,
        tokens: Iterable[Tuple[str, int, str]] = (),
        page: Optional[int] = None,
    ):  # (name, code, doc)
        self.array_alias: str = ""
        self.index = int(index)  # 1-based
        self.doc = doc
        self.tokens: Tuple[Tuple[str, int, str], ...] = tuple(tokens)
        self._attr_name: Optional[str] = None
        self.page: Optional[int] = int(page) if page is not None else None

    def __set_name__(self, owner, name):
        self._attr_name = name
        # Infer array alias from the owning class
        owner_alias = getattr(owner, "_array_alias", None)
        if not owner_alias:
            owner_alias = owner.__name__.lower()
        self.array_alias = str(owner_alias)

    # Helpers to read/write the owner's backing array (1-based indexing)
    def _get_code(self, instance: Any) -> float | int:
        return instance[self.index]

    def _set_code(self, instance: Any, value: float | int) -> None:
        instance[self.index] = value

    def __get__(self, instance, owner) -> Any:
        if instance is None:
            return self
        return BoundParam(
            owner_alias=self.array_alias.lower(),
            param_name=self._attr_name or "param",
            idx1=self.index,
            get_code=lambda: self._get_code(instance),
            set_code=lambda code: self._set_code(instance, code),
            doc=self.doc,
            tokens=self.tokens,
            page=self.page,
        )

    def __set__(self, instance, value):
        # Enumerated parameter: enforce allowed codes
        if self.tokens:
            if isinstance(value, BoundToken):
                code = int(value)
            elif isinstance(value, int):
                code = int(value)
            else:
                raise TypeError(
                    f"{self.array_alias}.{self._attr_name} expects a token or int, got {type(value).__name__}"
                )
            allowed = {code for _, code, _ in self.tokens}
            if code not in allowed:
                raise ValueError(
                    f"Illegal value for {self.array_alias}.{self._attr_name}: {code} not in {sorted(allowed)}"
                )
            self._set_code(instance, code)
            return
        # Non-enumerated: accept ints or floats, no validation
        if isinstance(value, (int, float)):
            instance[self.index] = value
            return
        raise TypeError(
            f"Unsupported assignment to {self.array_alias}({self.index}): {type(value)!r}"
        )


def param(*, index: int, page: Optional[int] = None):
    """Decorator to define a ParamDescriptor from an Enum-like nested class.

    Usage:
        @param(index=4)
        class message_level:
            "printing/verbosity level"
            errors_only = 1, "only error messages printed"
            ...
    """

    def wrapper(token_cls: type) -> ParamDescriptor:
        doc = (token_cls.__doc__ or "").strip()
        tokens: List[Tuple[str, int, str]] = []
        for name, val in token_cls.__dict__.items():
            if name.startswith("_"):
                continue
            if callable(val):
                continue
            code: int
            tdoc: str
            if isinstance(val, tuple) and len(val) >= 1:
                code = int(val[0])
                tdoc = str(val[1]) if len(val) >= 2 else ""
            elif isinstance(val, int):
                code = int(val)
                tdoc = ""
            else:
                # Skip unexpected shapes
                continue
            tokens.append((name, code, tdoc))
        return ParamDescriptor(index=index, doc=doc, tokens=tokens, page=page)

    return wrapper


class _OneBasedArray:
    """Simple 1-based indexable array wrapper for ints/floats."""

    def __init__(self, size: int, default: float | int = 0):
        self._data: List[float | int] = [default] * (size + 1)  # index 0 unused

    def __getitem__(self, idx1: int) -> float | int:
        if idx1 <= 0 or idx1 >= len(self._data):
            raise IndexError("1-based index out of range")
        return self._data[idx1]

    def __setitem__(self, idx1: int, value: float | int) -> None:
        if idx1 <= 0 or idx1 >= len(self._data):
            raise IndexError("1-based index out of range")
        self._data[idx1] = value


class ParamArray:
    """Base class for parameter arrays providing 1-based storage and aliasing."""

    _array_alias: Optional[str] = None
    _index_to_descriptor: Dict[int, ParamDescriptor] = {}

    def __init_subclass__(cls) -> None:
        if getattr(cls, "_array_alias", None) is None:
            cls._array_alias = cls.__name__.lower()
        # Build index->descriptor map for validation
        mapping: Dict[int, ParamDescriptor] = {}
        for name, attr in cls.__dict__.items():
            if isinstance(attr, ParamDescriptor):
                mapping[attr.index] = attr
        cls._index_to_descriptor = mapping

    def __init__(self, size: int = 64):
        self._array = _OneBasedArray(size, default=0)
        self._array_alias = type(self)._array_alias or type(self).__name__.lower()

    def __getitem__(self, idx1: int) -> float | int:
        return self._array[idx1]

    def __setitem__(self, idx1: int, value: float | int) -> None:
        desc = type(self)._index_to_descriptor.get(idx1)
        if desc and desc.tokens:
            # Enumerated param: enforce validation
            if isinstance(value, BoundToken):
                code = int(value)
            elif isinstance(value, int):
                code = int(value)
            else:
                raise TypeError(
                    f"{self._array_alias}.{desc._attr_name} expects a token or int, got {type(value).__name__}"
                )
            allowed = {c for _, c, _ in desc.tokens}
            if code not in allowed:
                raise ValueError(
                    f"Illegal value for {self._array_alias}.{desc._attr_name}: {code} not in {sorted(allowed)}"
                )
            self._array[idx1] = code
            return
        # Non-enumerated or unknown: store as-is (float or int)
        self._array[idx1] = value

    def __len__(self) -> int:
        # exclude index 0
        return max(0, len(self._array._data) - 1)

    def values(self) -> List[float | int]:
        return [self[i] for i in range(1, len(self) + 1)]


class RawArray(_OneBasedArray):
    """Unvalidated 1-based array (ints or floats depending on default)."""

    pass
