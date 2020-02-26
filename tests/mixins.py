import os


class FixtureMixin(object):
    @classmethod
    def setUpClass(cls):
        super(FixtureMixin, cls).setUpClass()
        test_dir = os.path.split(__file__)[0]
        cls.fixtures = os.path.join(test_dir, "fixtures")
