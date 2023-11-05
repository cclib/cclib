def str_contains_only(string, chars):
    """Checks if string contains only the specified characters.
    """
    return all([c in chars for c in string])
