# Load R libs ####
library(keras)
library(umap)
library(ggplot2)

# Load data ####
seu <- readRDS("../data/b_seu.rds")
expr <- seu@assays$RNA@scale.data
meta <- seu@meta.data

# Define train/test split ####
test <- sample(colnames(seu), round(ncol(seu)*0.1))
train <- setdiff(colnames(seu), test)

test <- rownames(meta)[meta$patient == "A"]
train <- setdiff(colnames(seu), test)

x_train <- expr[, train]
x_train <- t(x_train)

x_test <- expr[, test]
x_test <- t(x_test)

# Define outcome variable ####
outcome <- meta[, "ibrutinib_sensitivity"]
names(outcome) <- rownames(meta)
classes <- 0:4
names(classes) <- unique(outcome)
y_train <- to_categorical(classes[outcome[train]], 5)
y_test <- to_categorical(classes[outcome[test]], 5)

# Define model ####
input <- layer_input(shape = ncol(x_train))
embedding = input %>% 
  layer_dense(units = 64, activation = 'relu') %>% 
  layer_dropout(rate = 0.2) %>% 
  layer_dense(units = 16, activation = 'relu') %>%
  layer_dropout(rate = 0.1) %>% 
  layer_dense(units = 16, activation = 'relu') 
  
get_embedding = keras_model(input, embedding)

final_layer = embedding %>%
  layer_dense(units = 5, activation = 'softmax')

model = keras_model(input, final_layer)

model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = optimizer_adam(),
  metrics = c('accuracy')
)

# Train model ####
model %>% fit(
  x_train, y_train,
  batch_size = 64,
  epochs = 10,
  validation_split = 0.1
)

# Get predictions ####
tmp <- rbind(x_test, x_train)
predicted_classes <- predict(model, tmp)
predicted_classes <- apply(predicted_classes, 1, function(x) which(x == max(x)))
predicted_classes <- names(classes[predicted_classes])


table(predicted_classes[match(train, rownames(meta))], outcome[match(train, rownames(meta))])

table(predicted_classes[match(test, rownames(meta))], outcome[match(test, rownames(meta))])

embedding_coord <- predict(get_embedding, tmp)
uData <- umap(embedding_coord)

tmp <- rbind(data.frame(data = "test", meta[test, ]),
             data.frame(data = "train", meta[train,]))

aframe <- data.frame(uData$layout, tmp)
#aframe <- data.frame(embedding_coord, tmp)

ggplot(aframe, aes(X1, X2, color = ibrutinib_sensitivity)) +
  facet_wrap(~ data) +
  geom_point() +
  theme_bw()
