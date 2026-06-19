"""Custom exceptions for the ChemEase library."""


class ChemEaseError(Exception):
    """Base exception for all ChemEase errors."""

    pass


class InvalidMoleculeError(ChemEaseError):
    """Raised when an invalid or null molecule is encountered."""

    pass


class FetchError(ChemEaseError):
    """Raised when a molecule cannot be fetched from a database."""

    pass
