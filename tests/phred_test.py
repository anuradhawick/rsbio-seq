from rsbio_seq import ascii_to_phred, phred_to_ascii


def test_to_phred():
    # Test with typical input
    ascii_values = '!"#$%'
    expected_phred = [0, 1, 2, 3, 4]
    assert (
        ascii_to_phred(ascii_values) == expected_phred
    ), "Failed: ASCII to Phred conversion with typical input"

    # Test with empty string
    ascii_values_empty = ""
    expected_phred_empty = []
    assert (
        ascii_to_phred(ascii_values_empty) == expected_phred_empty
    ), "Failed: ASCII to Phred conversion with empty input"

    # Test with lower boundary edge case
    ascii_values_low = chr(33)
    expected_phred_low = [0]
    assert (
        ascii_to_phred(ascii_values_low) == expected_phred_low
    ), "Failed: ASCII to Phred conversion with lower boundary"


def test_to_ascii():
    # Test with typical input
    phred_scores = [0, 1, 2, 3, 4]
    expected_ascii = '!"#$%'
    assert (
        phred_to_ascii(phred_scores) == expected_ascii
    ), "Failed: Phred to ASCII conversion with typical input"

    # Test with empty list
    phred_scores_empty = []
    expected_ascii_empty = ""
    assert (
        phred_to_ascii(phred_scores_empty) == expected_ascii_empty
    ), "Failed: Phred to ASCII conversion with empty input"

    # Test with zero only
    phred_scores_zero = [0]
    expected_ascii_zero = "!"
    assert (
        phred_to_ascii(phred_scores_zero) == expected_ascii_zero
    ), "Failed: Phred to ASCII conversion with zero"
