# 更新 repo 的具体做法

## 本地 → GitHub (origin)

直接 `git push` 即可，GitHub 是 bare 仓库，commit + 文件一起更新。

```bash
git push origin devel-melt
```

## 本地 → mgt (SSH, 非 bare 仓库)

mgt 上的 `~/github/extrempy` 是**非 bare 仓库**（有 working tree），直接 `git push mgt` 只会更新 HEAD 指针和 object 数据库，**不会更新 working tree 和 index**。表现为：

- `git log` 看到新 commit ✓
- `git status` 仍显示旧文件为 "Changes to be committed" ✗
- 工作区文件内容仍为旧版本 ✗

### 正确做法（二选一）

**方式 A：先推 GitHub，再到 mgt 上 `git pull`**

```bash
# 本地
git push origin devel-melt

# mgt 上
cd ~/github/extrempy && git pull
```

**方式 B：本地同时推两个 remote，然后在 mgt 上 reset**

```bash
# 本地
git commit -m "..."
git push origin devel-melt
git push mgt devel-melt

# mgt 上
cd ~/github/extrempy && git reset --hard HEAD
```

推荐方式 A，更自然且安全。
