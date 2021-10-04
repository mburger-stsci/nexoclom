'''Exceptions used in nexoclom'''
class InputError(Exception):
    """Raised when a required parameter is not included"""
    def __init__(self, expression, message):
        self.expression = expression
        self.message = message

