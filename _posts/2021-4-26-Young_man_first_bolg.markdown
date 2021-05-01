---
layout: post
title: "年轻人的第一篇博客"
subtitle: 'First blog of a young man'
author: "louis"
mathjax: true
header-style: text
tags:
  - life
---

## 正文 

> 吃着火锅唱着歌, 就把博客搭好了

年轻人的第一篇博客,
hello world!

## 数学公式

行内公式: $g(x) = \frac{1 + x}{2x^2}$

行间公式: 
$$ f(x) = \frac{1}{2x} $$

- [ ] TODO: 研究mathjax的支持

## 代码块

```shell
echo hh -at 
```
```python
print("hello world")
for a in matrix:
  do_some_thing(a)

def fuc(fun, *args, **kwgs):
  def wrapper(fun, args, kwgs):
    print("I am upper wrapper")
    result = fun(args, kwgs)
    print("I am lower wrapper")
    return result
  return wrapper
    
@fuc
do_some_thing(matrix)
```

```java
public class A extends B implements C
{
  private String name;
  private int age;
  public A(String name, int age) {
    this.name = name;
    this.age  = age;
  }

  public void say() {
    System.out.println("I am " + name + "and I am " + age + " years old!")
  }

  public String toString() {
    return "<A: name: " + name + " age: " + age + ">"
  }

  public static void main(Stirng[] args) {
    A a = new A("a", 10)
    a.say()
  }
}
```