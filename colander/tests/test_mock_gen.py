from colander.mock_data_generation.utils import rand_nt

def test_rand_nt():
    for i in range(100):
        n = rand_nt()
        assert n in ["A", "C", "G", "T"]
