from textprogress import TextProgress
try:
    import qt
except ImportError:
    pass # import QtProgress will cause an error
else:
    from qtprogress import QtProgress
