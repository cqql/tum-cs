import Data.List

bitVectors 0 k = [[]]
bitVectors n 0 = [replicate n 0]
bitVectors n k
  | n == k = [replicate k 1]
  | otherwise =
      (map (0 :) (bitVectors (n - 1) k))
      ++ (map (1 :) (bitVectors (n - 1) (k - 1)))

showVector v = "[" ++ (intercalate ", " (map show v)) ++ "]"

vectors = bitVectors 15 9
sv = length vectors

rotate n xs = take (length xs) (drop n (cycle xs))

loneliness (a, b, c) = if b == 1 && not (a == 1 && c == 1) then 1 else 0

lonely v = sum (map loneliness (zip3 v (rotate 1 v) (rotate 2 v)))

p v = showVector v ++ " " ++ show (lonely v)

main = do
  putStrLn $ intercalate "\n" $ map p vectors
  putStrLn $ show ((toRational (sum (map lonely vectors))) / (toRational sv))
