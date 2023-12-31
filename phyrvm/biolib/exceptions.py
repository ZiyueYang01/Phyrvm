class PhyRvmException(Exception):
    def __init__(self, message=''):
        Exception.__init__(self, message)

class PhyRvmExit(Exception):
    def __init__(self, message=''):
        Exception.__init__(self, message)

class FileNotFound(PhyRvmException):
    def __init__(self, message=''):
        PhyRvmException.__init__(self, message)


class DirNotFound(PhyRvmException):
    def __init__(self, message=''):
        PhyRvmException.__init__(self, message)  
        
class ExternalException(PhyRvmException):
    def __init__(self, message=''):
        PhyRvmException.__init__(self, message)

