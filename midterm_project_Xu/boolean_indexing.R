# Test Dataframe
test <- data.frame(c(1,2,3,4,5,6), c(2,4,6,8,10,12))
colnames(test) <- c("one", "two")
# Making a mask that delete all data other than if test$one is 1
test_mask <- ifelse(test$one == 1, TRUE, FALSE)
# boolean indexing
test <- test[test_mask,]
