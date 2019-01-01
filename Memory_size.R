# Check your max allocated current v. memory size:

memory.size(max = TRUE)

# Save it:

ori_memory_size <- memory.size(max = TRUE)

# Increase it (carefully):

memory.size(ori_memory_size * 25)

# Increase factor gradually, until errors stop occurring. If computations still run slow, repeat procedure.

# From what I understand, the max v. memory size you assign to the R session is the size of your hardware RAM. You should be able to look that up in the properties of your system/computer.