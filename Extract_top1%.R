table = read.csv("migrant_vs_resident_Final_stats.csv")
test1 = table[order(-table$West_vs_resident_Fst), ]
westvresident = test1[1:236, ]
View(westvresident)
test2 = table[order(-table$East_vs_resident_Fst), ]
eastvresident = test2[1:236, ]
test3 = table[order(-table$migrant_vs_resident_Fst), ]
migrantvresident = test3[1:236, ]

write.csv(westvresident, "West_vs_Resident.csv")
write.csv(eastvresident, "East_vs_Resident.csv")
write.csv(migrantvresident, "migrant_vs_Resident.csv")