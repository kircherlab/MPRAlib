class MPRAlibException(Exception):
    """
    MPRAlib error class for specific exceptions.

    Attributes:
        message (str): A description of the error.
    """

    def __init__(self, message: str):
        super().__init__(message)
        self.message = message

    def __str__(self):
        return f"{self.__class__.__name__}: {self.message}"


class IOException(MPRAlibException):
    """
    Exception raised for IO-related errors.

    Attributes:
        message (str): A description of the IO error.
    """

    def __init__(self, message: str):
        super().__init__(message)


class SequenceDesignException(IOException):
    """
    Exception raised for errors related to sequence design file.

    Attributes:
        message (str): A description of the sequence design error.
    """

    def __init__(self, column: str, file_path: str):
        super().__init__(f"Column {column} has a wrong format in file {file_path}.")
