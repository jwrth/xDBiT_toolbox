class ModuleNotFoundOnWindows(ModuleNotFoundError):
    '''
    Code from https://github.com/theislab/scib/blob/main/scib/exceptions.py
    Information about structure: https://careerkarma.com/blog/python-super/
    '''

    def __init__(self, exception):
        self.message = f"\n{exception.name} is not installed. " \
                       "This package could be problematic to install on Windows."
        super().__init__(self.message)

