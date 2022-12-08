from unittest import TestCase

from mtbc import MtbcRandom


class TestMtbcRandom(TestCase):

    def test_construct_search_request(self):
        mtbc = MtbcRandom()
        mtbc_78331 = MtbcRandom(select_mycobacterium_mungi=True)
        mtbc_78331_1844474 = MtbcRandom(select_mycobacterium_mungi=True, select_mycobacterium_canettii=True)
        self.assertEqual(mtbc.req, 'txid77643[ORGN]')
        self.assertEqual(mtbc_78331.req, 'txid1844474[ORGN]')
        self.assertEqual(mtbc_78331_1844474.req, 'txid78331[ORGN] OR txid1844474[ORGN]')
    #
    # def test_mtbc_random_search_id(self):
    #     self.fail()
    #
    # def test_mtbc_random_get_accession_number(self):
    #     self.fail()
