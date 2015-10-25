#include <stdlib.h>
#include <stdio.h>
#include <time.h>

typedef struct Node {
  struct Node* left;
  struct Node* right;
  int nodeValue;
} Node;

Node* makeRandomTree(int numberOfNodes) {
  if (numberOfNodes == 0) {
    return NULL;
  } else {
    Node* tree = malloc(sizeof(Node));
    tree->nodeValue = rand() % 100;

    if (numberOfNodes > 1) {
      int n = numberOfNodes - 1;
      int leftNum = rand() % n;
      int rightNum = n - leftNum;

      tree->left = makeRandomTree(leftNum);
      tree->right = makeRandomTree(rightNum);
    }

    return tree;
  }
}

int main() {
  srand(time(NULL));

  int n = rand();

  printf("Generate tree with n = %d\n", n);

  makeRandomTree(n);

  return 0;
}
