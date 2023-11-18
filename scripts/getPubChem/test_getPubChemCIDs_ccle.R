source("newfunctions.R")
# # Test case 1: simple test
name1 <- "AcetyLchOline"
expected_output1 <- "ACETYLCHOLINE"
assertthat::assert_that(cleanCharacterStrings(name1) == expected_output1) 

# Test case 2: Test removing round brackets and contents
name2 <- "Acetylcholine (ACh)"
expected_output2 <- "ACETYLCHOLINE"
assertthat::assert_that(cleanCharacterStrings(name2) == expected_output2)

# Test case 3: Test removing square brackets and contents
name3 <- "Acetylcholine [ACh]"
expected_output3 <- "ACETYLCHOLINE"
assertthat::assert_that(cleanCharacterStrings(name3) == expected_output3)

# Test case 4: Test removing curly brackets and contents
name4 <- "Acetylcholine {ACh}"
expected_output4 <- "ACETYLCHOLINE"
assertthat::assert_that(cleanCharacterStrings(name4) == expected_output4)

# Test case 5: Test converting entire string to uppercase
name5 <- "Acetylcholine"
expected_output5 <- "ACETYLCHOLINE"
assertthat::assert_that(cleanCharacterStrings(name5) == expected_output5)# Test case 1: Test removing micro symbol

# Test case 6: Test removing comma
name6 <- "Acetylcholine, chloride"
expected_output6 <- "ACETYLCHOLINECHLORIDE"
assertthat::assert_that(cleanCharacterStrings(name6) == expected_output6)

# Test case 7: Test removing semicolon
name7 <- "Acetylcholine; chloride"
expected_output7 <- "ACETYLCHOLINECHLORIDE"
assertthat::assert_that(cleanCharacterStrings(name7) == expected_output7)

# Test case 8: Test removing colon
name8 <- "Acetylcholine: chloride"
expected_output8 <- "ACETYLCHOLINECHLORIDE"
assertthat::assert_that(cleanCharacterStrings(name8) == expected_output8)

# Test case 9: Test removing hyphen
name9 <- "Acetylcholine - chloride"
expected_output9 <- "ACETYLCHOLINECHLORIDE"
assertthat::assert_that(cleanCharacterStrings(name9) == expected_output9)

# Test case 10: Test removing plus sign
name10 <- "Acetylcholine + chloride"
expected_output10 <- "ACETYLCHOLINECHLORIDE"
assertthat::assert_that(cleanCharacterStrings(name10) == expected_output10)

# Test case 11: Test removing asterisk
name11 <- "Acetylcholine * chloride"
expected_output11 <- "ACETYLCHOLINECHLORIDE"
assertthat::assert_that(cleanCharacterStrings(name11) == expected_output11)

# Test case 12: Test removing percent sign
name12 <- "Acetylcholine % chloride"
expected_output12 <- "ACETYLCHOLINECHLORIDE"
assertthat::assert_that(cleanCharacterStrings(name12) == expected_output12)

# Test case 13: Test removing dollar sign
name13 <- "Acetylcholine $ chloride"
expected_output13 <- "ACETYLCHOLINECHLORIDE"
assertthat::assert_that(cleanCharacterStrings(name13) == expected_output13)

# Test case 14: Test removing hash sign
name14 <- "Acetylcholine # chloride"
expected_output14 <- "ACETYLCHOLINECHLORIDE"
assertthat::assert_that(cleanCharacterStrings(name14) == expected_output14)

# Test case 15: Test removing underscore
name15 <- "Acetylcholine_chloride"
expected_output15 <- "ACETYLCHOLINECHLORIDE"
assertthat::assert_that(cleanCharacterStrings(name15) == expected_output15)

# Test case 16: Test removing space
name16 <- "Acetylcholine chloride"
expected_output16 <- "ACETYLCHOLINECHLORIDE"
assertthat::assert_that(cleanCharacterStrings(name16) == expected_output16)

# Test case 17: Test removing round brackets and contents with spaces
name17 <- "Acetylcholine (ACh) chloride"
expected_output17 <- "ACETYLCHOLINECHLORIDE"
assertthat::assert_that(cleanCharacterStrings(name17) == expected_output17)

# Test case 18: Test removing square brackets and contents with spaces
name18 <- "Acetylcholine [ACh] chloride"
expected_output18 <- "ACETYLCHOLINECHLORIDE"
assertthat::assert_that(cleanCharacterStrings(name18) == expected_output18)

# Test case 19: Test removing curly brackets and contents with spaces
name19 <- "Acetylcholine {ACh} chloride"
expected_output19 <- "ACETYLCHOLINECHLORIDE"
assertthat::assert_that(cleanCharacterStrings(name19) == expected_output19)

# Test case 20: Test removing invalid character
name20 <- "Acetylcholine @ chloride"
expected_output20 <- "ACETYLCHOLINE@CHLORIDE"
assertthat::assert_that(cleanCharacterStrings(name20) == expected_output20)